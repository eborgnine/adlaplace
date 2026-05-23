#' Parameter-block row counts for an \code{ad_model}
#'
#' @param x An \code{ad_model} object.
#' @return Named list with \code{n_beta}, \code{n_gamma}, and \code{n_theta}.
#' @export
setGeneric("sizes", function(x, ...) standardGeneric("sizes"))

#' @rdname sizes
#' @export
setMethod("sizes", "ad_model", function(x) {
  list(
    n_beta = nrow(x@beta_map),
    n_gamma = nrow(x@gamma_map),
    n_theta = nrow(x@theta_map)
  )
})

#' @keywords internal
validate_ad_model_dims <- function(y, A, X) {
  n_obs <- length(y)
  if (n_obs == 0L) {
    if (!is.null(A) && (nrow(A) != 0L || ncol(A) != 0L)) {
      stop("y is empty but A is non-empty")
    }
    if (!is.null(X) && (nrow(X) != 0L || ncol(X) != 0L)) {
      stop("y is empty but X is non-empty")
    }
    return(invisible(NULL))
  }
  if (!is.null(A) && nrow(A) != n_obs) {
    stop("nrow(A) (", nrow(A), ") must match length(y) (", n_obs, ")")
  }
  if (!is.null(X) && nrow(X) != n_obs) {
    stop("nrow(X) (", nrow(X), ") must match length(y) (", n_obs, ")")
  }
  invisible(NULL)
}

#' @keywords internal
validate_ad_model_maps <- function(model, kind) {
  layout <- sizes(model)
  if (identical(kind, "parameters")) {
    if (layout$n_theta == 0L || ncol(model@theta_map) == 0L) {
      stop("parameters shards require a non-empty theta_map")
    }
  } else if (identical(kind, "random")) {
    if (ncol(model@beta_map) != 0L) {
      stop("random shards require beta_map with zero columns")
    }
    if (layout$n_gamma == 0L) {
      stop("random shards require a non-empty gamma_map")
    }
  } else if (identical(kind, "observations")) {
    if (length(model@y) == 0L) {
      stop("observation shards require non-empty y on ad_model")
    }
  }
  invisible(layout)
}

#' @keywords internal
infer_n_total <- function(info = list(), config = NULL) {
  if (!is.null(info$parameters) && !is.null(info$gamma)) {
    return(as.integer(nrow(info$parameters) + nrow(info$gamma)))
  }
  if (!is.null(config)) {
    n_beta <- if (!is.null(config$beta)) length(config$beta) else NA_integer_
    n_gamma <- if (!is.null(config$gamma)) length(config$gamma) else NA_integer_
    n_theta <- if (!is.null(config$theta)) length(config$theta) else NA_integer_
    if (!any(is.na(c(n_beta, n_gamma, n_theta))) &&
        n_beta + n_gamma + n_theta > 0L) {
      return(as.integer(n_beta + n_gamma + n_theta))
    }
  }
  NA_integer_
}

#' @keywords internal
infer_n_theta <- function(n_beta,
                          n_gamma,
                          theta_map = NULL,
                          info = list(),
                          config = NULL) {
  if (!is.null(theta_map)) {
    return(as.integer(nrow(theta_map)))
  }
  if (!is.null(info$theta)) {
    return(as.integer(nrow(info$theta)))
  }
  n_total <- infer_n_total(info = info, config = config)
  if (!is.na(n_total) && !is.na(n_gamma) && !is.na(n_beta)) {
    n_theta <- as.integer(n_total - n_beta - n_gamma)
    if (n_theta < 0L) {
      stop(
        "cannot infer n_theta: n_total (", n_total,
        ") < n_beta (", n_beta, ") + n_gamma (", n_gamma, ")",
        call. = FALSE
      )
    }
    return(n_theta)
  }
  stop(
    "cannot infer n_theta: provide theta_map, info$theta, or enough ",
    "information to form n_total (info$parameters + info$gamma, or ",
    "config$beta + config$gamma + config$theta) together with n_beta ",
    "and n_gamma",
    call. = FALSE
  )
}

#' @keywords internal
validate_config_layout <- function(model, config, kind = NULL) {
  layout <- sizes(model)
  if (length(config$beta) != layout$n_beta) {
    stop(
      "length(config$beta) (", length(config$beta),
      ") must match nrow(beta_map) (", layout$n_beta, ")"
    )
  }
  if (is.null(kind) || !identical(kind, "random")) {
    if (length(config$gamma) != layout$n_gamma) {
      stop(
        "length(config$gamma) (", length(config$gamma),
        ") must match nrow(gamma_map) (", layout$n_gamma, ")"
      )
    }
  }
  n_theta_config <- if (!is.null(config$theta)) length(config$theta) else NA_integer_
  if (is.na(n_theta_config) || n_theta_config == 0L) {
    n_theta_config <- infer_n_theta(
      n_beta = layout$n_beta,
      n_gamma = layout$n_gamma,
      info = model@info,
      config = config
    )
  }
  if (n_theta_config != layout$n_theta) {
    stop(
      "length(config$theta) (", n_theta_config,
      ") must match nrow(theta_map) (", layout$n_theta, ")"
    )
  }
  invisible(layout)
}

#' Coerce a sparse matrix to \code{dgCMatrix}.
#' @keywords internal
as_dgC <- function(m) {
  methods::as(methods::as(m, "CsparseMatrix"), "dgCMatrix")
}

#' Coerce a sparse matrix to \code{ngCMatrix} (Matrix >= 1.6).
#' @keywords internal
as_ngC <- function(m) {
  if (inherits(m, "ngCMatrix")) {
    return(m)
  }
  methods::as(methods::as(as_dgC(m), "nMatrix"), "ngCMatrix")
}

#' Coerce map shorthand: length-1 \code{nrow} (zero columns) or length-2
#' \code{c(row, nrow)} (one column, structural 1 at row).
#' @keywords internal
coerce_map_shorthand <- function(x, what) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.atomic(x)) {
    return(x)
  }
  if (length(x) == 1L) {
    n <- as.integer(x)
    if (is.na(n)) {
      return(x)
    }
    if (n < 0L) {
      stop(what, ": nrow must be non-negative", call. = FALSE)
    }
    return(Matrix::Matrix(nrow = n, ncol = 0L))
  }
  if (length(x) == 2L) {
    row <- as.integer(x[1])
    n <- as.integer(x[2])
    if (is.na(row) || is.na(n)) {
      stop(what, ": shorthand c(row, nrow) must be finite integers", call. = FALSE)
    }
    if (n < 0L) {
      stop(what, ": nrow must be non-negative", call. = FALSE)
    }
    if (n == 0L) {
      return(Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        dims = c(0L, 1L),
        giveCsparse = TRUE
      ))
    }
    if (row < 1L || row > n) {
      stop(
        what, ": row index (", row, ") must be between 1 and nrow (", n, ")",
        call. = FALSE
      )
    }
    return(Matrix::sparseMatrix(
      i = row,
      j = 1L,
      x = 1,
      dims = c(n, 1L),
      giveCsparse = TRUE
    ))
  }
  x
}

#' @keywords internal
validate_map_at_most_one_per_row_col <- function(m, what) {
  if (nrow(m) == 0L || ncol(m) == 0L) {
    return(invisible(NULL))
  }
  rs <- Matrix::rowSums(m)
  cs <- Matrix::colSums(m)
  if (any(rs > 1L)) {
    stop(what, ": more than one nonzero per row", call. = FALSE)
  }
  if (any(cs > 1L)) {
    stop(what, ": more than one nonzero per column", call. = FALSE)
  }
  invisible(NULL)
}

#' Transpose of a design matrix for \code{ATp} / \code{XTp} slots.
#' @keywords internal
design_Tp <- function(M, n_obs, what = "M") {
  if (is.null(M)) {
    return(as_dgC(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(0L, n_obs),
      giveCsparse = TRUE
    )))
  }
  if (nrow(M) != n_obs) {
    stop(
      "nrow(", what, ") (", nrow(M), ") must match length(y) (", n_obs, ")",
      call. = FALSE
    )
  }
  as_dgC(Matrix::t(M))
}

#' Build \code{ad_model} from designs and \code{config} layout sizes.
#'
#' Used by examples with a plain data list (\code{y}, \code{A}, \code{X}).
#'
#' @keywords internal
ad_model_from_config_matrices <- function(y, A, X, config, theta_local_row = 0L) {
  n_gamma <- length(config$gamma)
  n_theta <- length(config$theta)
  gamma_map <- if (n_gamma > 0L) {
    Matrix::Matrix(nrow = n_gamma, ncol = 0L)
  } else {
    NULL
  }
  theta_map <- if (n_theta > 0L) {
    Matrix::sparseMatrix(
      i = as.integer(theta_local_row),
      j = 0L,
      x = 1,
      dims = c(n_theta, 1L),
      index1 = FALSE,
      giveCsparse = TRUE
    )
  } else {
    NULL
  }
  ad_model(
    y = y,
    A = A,
    X = X,
    gamma_map = gamma_map,
    theta_map = theta_map
  )
}

#' @describeIn model_setup Construct an ad_model object
#'
#' @param y Response vector (default empty).
#' @param A Random-effects design matrix (\code{nrow(A) = length(y)}; any
#'   \code{Matrix} class).
#' @param X Fixed-effects design matrix (\code{nrow(X) = length(y)}; any
#'   \code{Matrix} class).
#' @param beta_map Beta parameter map (any \code{Matrix}; coerced to
#'   \code{ngCMatrix}), or a shorthand: length-1 integer \code{n} gives
#'   \code{Matrix::Matrix(nrow = n, ncol = 0)}; length-2 \code{c(row, n)} gives
#'   one column with a structural \code{1} at row \code{row} (1-based). If
#'   \code{NULL}, defaults to \code{Matrix::Diagonal(ncol(X))} when
#'   \code{ncol(X) > 0}, otherwise an empty matrix.
#' @param gamma_map Gamma parameter map (any \code{Matrix}; coerced to
#'   \code{ngCMatrix}; at most one nonzero per row and column), or the same
#'   length-1 / length-2 shorthands as \code{beta_map}. If \code{NULL}, defaults
#'   to \code{Matrix::Matrix(nrow = ncol(A), ncol = 0)} when \code{ncol(A) > 0},
#'   otherwise an empty matrix.
#' @param theta_map Theta parameter map (any \code{Matrix}; coerced to
#'   \code{ngCMatrix}), or the same length-1 / length-2 shorthands as
#'   \code{beta_map}. If \code{NULL}, \code{n_theta} is inferred via
#'   \code{infer_n_theta()} from \code{info} and map dimensions, then defaults
#'   to \code{Matrix::Diagonal(n_theta)}.
#' @param elgmMatrix Optional ELGM map (any \code{Matrix}; coerced to
#'   \code{ngCMatrix}). If \code{NULL}, defaults to \code{Matrix::Matrix(nrow = length(y), ncol = 0)}.
#' @param data Original data frame (default empty).
#' @param terms Model terms (default empty list).
#' @param info Parameter info (default empty list).
#' @export
ad_model <- function(y = numeric(0),
                     A = NULL,
                     X = NULL,
                     beta_map = NULL,
                     gamma_map = NULL,
                     theta_map = NULL,
                     elgmMatrix = NULL,
                     data = NULL,
                     terms = list(),
                     info = list()) {
  n_obs <- length(y)
  validate_ad_model_dims(y, A, X)
  ATp <- design_Tp(A, n_obs, "A")
  XTp <- design_Tp(X, n_obs, "X")
  if (is.null(data)) {
    data <- data.frame()
  }
  beta_map <- coerce_map_shorthand(beta_map, "beta_map")
  gamma_map <- coerce_map_shorthand(gamma_map, "gamma_map")
  theta_map <- coerce_map_shorthand(theta_map, "theta_map")
  if (is.null(beta_map)) {
    n_beta <- if (!is.null(X) && ncol(X) > 0L) ncol(X) else 0L
    beta_map <- if (n_beta > 0L) {
      Matrix::Diagonal(n_beta)
    } else {
      Matrix::Matrix(nrow = 0L, ncol = 0L)
    }
  }
  if (is.null(gamma_map)) {
    n_gamma <- if (!is.null(A) && ncol(A) > 0L) ncol(A) else 0L
    gamma_map <- if (n_gamma > 0L) {
      Matrix::Matrix(nrow = n_gamma, ncol = 0L)
    } else {
      Matrix::Matrix(nrow = 0L, ncol = 0L)
    }
  }
  if (is.null(theta_map)) {
    n_theta <- infer_n_theta(
      n_beta = nrow(beta_map),
      n_gamma = nrow(gamma_map),
      info = info
    )
    theta_map <- if (n_theta > 0L) {
      Matrix::Diagonal(n_theta)
    } else {
      Matrix::Matrix(nrow = 0L, ncol = 0L)
    }
  }
  if (is.null(elgmMatrix)) {
    elgmMatrix <- Matrix::Matrix(nrow = n_obs, ncol = 0L)
  }
  beta_map <- as_ngC(beta_map)
  gamma_map <- as_ngC(gamma_map)
  theta_map <- as_ngC(theta_map)
  elgmMatrix <- as_ngC(elgmMatrix)
  validate_map_at_most_one_per_row_col(beta_map, "beta_map")
  validate_map_at_most_one_per_row_col(gamma_map, "gamma_map")
  methods::new(
    "ad_model",
    data = data,
    y = as.numeric(y),
    ATp = ATp,
    XTp = XTp,
    beta_map = beta_map,
    gamma_map = gamma_map,
    theta_map = theta_map,
    elgmMatrix = elgmMatrix,
    terms = terms,
    info = info
  )
}
