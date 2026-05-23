#' Model object from \code{model_setup()} or \code{ad_model()}
#'
#' @slot data Original model \code{data.frame}.
#' @slot y Response vector.
#' @slot ATp Transpose of random-effects design matrix.
#' @slot XTp Transpose of fixed-effects design matrix.
#' @slot beta_map Beta parameter map (\code{ngCMatrix}).
#' @slot gamma_map Gamma parameter map (\code{ngCMatrix}).
#' @slot theta_map Theta parameter map (\code{ngCMatrix}).
#' @slot elgmMatrix Optional exposure–lag map (\code{ngCMatrix}; empty by default).
#' @slot terms Named list of model term objects.
#' @slot info List with \code{beta}, \code{gamma}, \code{theta}, \code{parameters}.
#' @exportClass ad_model
setClass(
  "ad_model",
  slots = c(
    data = "data.frame",
    y = "numeric",
    ATp = "Matrix",
    XTp = "Matrix",
    beta_map = "ngCMatrix",
    gamma_map = "ngCMatrix",
    theta_map = "ngCMatrix",
    elgmMatrix = "ngCMatrix",
    terms = "list",
    info = "list"
  )
)

#' @keywords internal
empty_ngC <- function() {
  methods::as(
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(0L, 0L),
      giveCsparse = TRUE
    ),
    "ngCMatrix"
  )
}

#' @keywords internal
default_elgm <- function(n_obs) {
  methods::as(
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(as.integer(n_obs), 0L),
      giveCsparse = TRUE
    ),
    "ngCMatrix"
  )
}

#' @keywords internal
empty_matrix <- function(nrow = 0L, ncol = 0L) {
  methods::as(
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(as.integer(nrow), as.integer(ncol)),
      giveCsparse = TRUE
    ),
    "dMatrix"
  )
}

#' @keywords internal
beta_map_zero_cols <- function(n_beta) {
  if (n_beta <= 0L) {
    return(empty_ngC())
  }
  methods::as(
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(as.integer(n_beta), 0L),
      giveCsparse = TRUE
    ),
    "ngCMatrix"
  )
}

#' @keywords internal
gamma_map_zero_cols <- function(n_gamma) {
  if (n_gamma <= 0L) {
    return(empty_ngC())
  }
  methods::as(
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(as.integer(n_gamma), 0L),
      giveCsparse = TRUE
    ),
    "ngCMatrix"
  )
}

#' Layout dimensions from an ad_model
#'
#' @param model An \code{ad_model} object.
#' @return List with \code{n_beta}, \code{n_gamma}, \code{n_theta}.
#' @keywords internal
ad_model_layout <- function(model) {
  if (!is(model, "ad_model")) {
    stop("model must be an ad_model object")
  }
  list(
    n_beta = nrow(model@beta_map),
    n_gamma = nrow(model@gamma_map),
    n_theta = nrow(model@theta_map)
  )
}

#' @keywords internal
validate_ad_model_dims <- function(y, ATp, XTp) {
  n_obs <- length(y)
  if (n_obs == 0L) {
    if (ncol(ATp) != 0L || ncol(XTp) != 0L) {
      stop("y is empty but ATp or XTp has columns")
    }
    return(invisible(NULL))
  }
  if (ncol(ATp) != n_obs) {
    stop("ncol(ATp) (", ncol(ATp), ") must match length(y) (", n_obs, ")")
  }
  if (ncol(XTp) != n_obs) {
    stop("ncol(XTp) (", ncol(XTp), ") must match length(y) (", n_obs, ")")
  }
  invisible(NULL)
}

#' @keywords internal
validate_ad_model_maps <- function(model, kind) {
  layout <- ad_model_layout(model)
  if (identical(kind, "parameters")) {
    if (ncol(model@beta_map) != 0L) {
      stop("parameters shards require beta_map with zero columns")
    }
    if (ncol(model@gamma_map) != 0L) {
      stop("parameters shards require gamma_map with zero columns")
    }
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
validate_config_layout <- function(model, config, kind = NULL) {
  layout <- ad_model_layout(model)
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
  if (length(config$theta) != layout$n_theta) {
    stop(
      "length(config$theta) (", length(config$theta),
      ") must match nrow(theta_map) (", layout$n_theta, ")"
    )
  }
  invisible(layout)
}

#' Coerce a sparse matrix to \code{ngCMatrix} (Matrix >= 1.6).
#' @keywords internal
as_ngC <- function(m) {
  cm <- methods::as(m, "CsparseMatrix")
  methods::as(methods::as(cm, "nMatrix"), "ngCMatrix")
}

#' Build \code{ad_model} from designs and \code{config} layout sizes.
#'
#' Used by examples with a plain data list (\code{y}, \code{ATp}, \code{XTp}).
#'
#' @keywords internal
ad_model_from_config_matrices <- function(y, ATp, XTp, config, theta_local_row = 0L) {
  n_beta <- length(config$beta)
  n_gamma <- length(config$gamma)
  n_theta <- length(config$theta)
  if (n_beta > 0L) {
    beta_map <- as_ngC(
      Matrix::sparseMatrix(
        i = seq.int(0L, n_beta - 1L),
        j = seq.int(0L, n_beta - 1L),
        dims = c(n_beta, n_beta),
        index1 = FALSE,
        giveCsparse = TRUE
      )
    )
  } else {
    beta_map <- empty_ngC()
  }
  if (n_gamma > 0L) {
    gamma_map <- as_ngC(
      Matrix::sparseMatrix(
        i = integer(0),
        j = integer(0),
        dims = c(n_gamma, 0L),
        giveCsparse = TRUE
      )
    )
  } else {
    gamma_map <- empty_ngC()
  }
  if (n_theta > 0L) {
    theta_map <- as_ngC(
      Matrix::sparseMatrix(
        i = as.integer(theta_local_row),
        j = 0L,
        x = 1,
        dims = c(n_theta, 1L),
        index1 = FALSE,
        giveCsparse = TRUE
      )
    )
  } else {
    theta_map <- empty_ngC()
  }
  ad_model(
    y = y,
    ATp = ATp,
    XTp = XTp,
    beta_map = beta_map,
    gamma_map = gamma_map,
    theta_map = theta_map
  )
}

#' @describeIn model_setup Construct an ad_model object
#'
#' @param y Response vector (default empty).
#' @param ATp Transpose of random-effects design matrix (default empty).
#' @param XTp Transpose of fixed-effects design matrix (default empty).
#' @param beta_map Beta parameter map (\code{ngCMatrix}).
#' @param gamma_map Gamma parameter map (\code{ngCMatrix}).
#' @param theta_map Theta parameter map (\code{ngCMatrix}).
#' @param elgmMatrix Optional ELGM map (\code{ngCMatrix}).
#' @param data Original data frame (default empty).
#' @param terms Model terms (default empty list).
#' @param info Parameter info (default empty list).
#' @export
ad_model <- function(y = numeric(0),
                     ATp = NULL,
                     XTp = NULL,
                     beta_map = empty_ngC(),
                     gamma_map = empty_ngC(),
                     theta_map = empty_ngC(),
                     elgmMatrix = NULL,
                     data = NULL,
                     terms = list(),
                     info = list()) {
  if (is.null(ATp)) {
    ATp <- empty_matrix(length(y), 0L)
  }
  if (is.null(XTp)) {
    XTp <- empty_matrix(length(y), 0L)
  }
  validate_ad_model_dims(y, ATp, XTp)
  n_obs <- length(y)
  if (is.null(elgmMatrix)) {
    elgmMatrix <- default_elgm(n_obs)
  }
  if (is.null(data)) {
    data <- data.frame()
  }
  beta_map <- methods::as(beta_map, "ngCMatrix")
  gamma_map <- methods::as(gamma_map, "ngCMatrix")
  theta_map <- methods::as(theta_map, "ngCMatrix")
  elgmMatrix <- methods::as(elgmMatrix, "ngCMatrix")
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

#' Maps-only ad_model for a random-effect shard
#'
#' @keywords internal
ad_model_random_maps <- function(gamma_map,
                                   theta_map,
                                   n_beta = 0L,
                                   n_theta = NULL) {
  if (is.null(n_theta)) {
    n_theta <- nrow(theta_map)
  }
  ad_model(
    beta_map = beta_map_zero_cols(n_beta),
    gamma_map = gamma_map,
    theta_map = theta_map
  )
}

#' Parameters-shard view of an ad_model (zero-column beta/gamma maps)
#'
#' @keywords internal
ad_model_parameters_view <- function(model) {
  layout <- ad_model_layout(model)
  ad_model(
    y = model@y,
    ATp = model@ATp,
    XTp = model@XTp,
    beta_map = beta_map_zero_cols(layout$n_beta),
    gamma_map = gamma_map_zero_cols(layout$n_gamma),
    theta_map = model@theta_map,
    elgmMatrix = model@elgmMatrix,
    data = model@data,
    terms = model@terms,
    info = model@info
  )
}
