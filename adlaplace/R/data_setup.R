#' Apply log transform to selected theta columns
#'
#' When \code{active} is \code{TRUE}, logs values in \code{cols} for rows where
#' \code{theta_info$transform} is \code{TRUE}. Missing or \code{NA} \code{transform}
#' entries default to \code{TRUE}.
#'
#' @param theta_info Data frame from \code{info$theta} (see \code{\link{data_setup}}).
#' @param cols Character vector of column names to transform (e.g. \code{"init"}).
#' @param active Logical; when \code{FALSE}, returns \code{theta_info} unchanged.
#' @return A copy of \code{theta_info} with selected entries logged.
#' @export
apply_theta_log <- function(theta_info, cols = "init", active = TRUE) {
  if (!active || nrow(theta_info) == 0L) {
    return(theta_info)
  }
  out <- theta_info
  transform <- if ("transform" %in% names(out)) {
    out$transform
  } else {
    rep(TRUE, nrow(out))
  }
  transform[is.na(transform)] <- TRUE
  idx <- which(transform)
  if (length(idx) == 0L) {
    return(out)
  }
  for (col in cols) {
    if (col %in% names(out)) {
      out[idx, col] <- log(out[idx, col])
    }
  }
  out
}

#' @keywords internal
ensure_theta_transform <- function(theta_setup) {
  if (nrow(theta_setup) == 0L) {
    theta_setup$transform <- logical(0)
    return(theta_setup)
  }
  if (!"transform" %in% names(theta_setup)) {
    theta_setup$transform <- TRUE
  } else {
    theta_setup$transform[is.na(theta_setup$transform)] <- TRUE
  }
  theta_setup
}

#' Build design matrices and parameter metadata from a formula
#'
#' @param formula Model formula or named list of term objects.
#' @param data Data frame.
#' @param verbose Print parsing messages.
#' @return List with \code{terms}, \code{y}, \code{A}, \code{X}, \code{data}, \code{info}.
#'   The \code{info$theta} data frame includes a logical \code{transform} column:
#'   \code{TRUE} means the parameter is stored on the log scale when
#'   \code{config$transform_theta} is \code{TRUE} (default for all rows unless a
#'   model term sets otherwise). \code{info$parameters} stacks \code{beta} and
#'   \code{theta} rows with the same \code{transform} column (\code{FALSE} for
#'   fixed effects); \code{init}, \code{lower}, and \code{upper} are already on
#'   the optimization (log) scale where \code{transform} is \code{TRUE}.
#' @export
data_setup <- function(formula, data, verbose = FALSE) {
  formula_in <- formula
  if (methods::is(formula, "formula")) {
    terms <- collect_terms(formula, verbose = verbose)
  } else {
    terms <- formula
    if (inherits(terms, "model")) {
      terms <- list(terms)
    }
    inherits_seq <- unlist(lapply(terms, inherits, what = "model"))
    if (!all(inherits_seq)) {
      warning(
        "formula must be of class formula or a list of objects which inherit class model"
      )
    }
  }

  terms <- lapply(terms, add_by_levels, data)

  design_list <- lapply(terms, design, data = data)

  terms_with_gamma <- sapply(terms, methods::slot, "type") == "random"
  terms_with_beta <- sapply(terms, methods::slot, "type") == "fixed"

  if (any(terms_with_beta)) {
    x_matrix <- do.call(cbind, design_list[terms_with_beta])
  } else {
    x_matrix <- matrix(nrow = nrow(data), ncol = 0)
  }

  if (any(terms_with_gamma)) {
    a_matrix <- do.call(cbind, design_list[terms_with_gamma])
  } else {
    a_matrix <- matrix(nrow = nrow(data), ncol = 0)
  }

  theta_info_list <- Filter(Negate(is.null), lapply(terms, theta_info))
  theta_info_list <- lapply(theta_info_list, function(df) {
    if (!"transform" %in% names(df)) {
      df$transform <- TRUE
    } else {
      df$transform[is.na(df$transform)] <- TRUE
    }
    df
  })
  if (length(theta_info_list) == 0L) {
    theta_setup <- data.frame(
      term = character(),
      model = character(),
      label = character(),
      init = numeric(),
      lower = numeric(),
      upper = numeric(),
      parscale = numeric(),
      type = factor(levels = .type_factor_levels),
      transform = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    theta_setup <- do.call(rbind, theta_info_list)
  }
  theta_setup <- ensure_theta_transform(theta_setup)

  theta_setup$id <- seq.int(0, length.out = nrow(theta_setup))

  beta_setup <- do.call(rbind, lapply(terms, beta_info, data = data))

  random_info_list <- Filter(Negate(is.null), lapply(terms, random_info, data = data))
  if (length(random_info_list) == 0L) {
    gamma_setup <- data.frame(
      term = character(),
      model = character(),
      label = character(),
      gamma_label = character(),
      stringsAsFactors = FALSE
    )
  } else {
    gamma_setup <- do.call(rbind, random_info_list)
  }
  gamma_setup$id <- seq.int(0, length.out = nrow(gamma_setup))
  if (nrow(gamma_setup) > 0L && nrow(theta_setup) > 0L) {
    gamma_setup$theta_id <- theta_setup[match(
      gamma_setup$label, theta_setup$label
    ), "id"]
  }

  if (ncol(x_matrix) > 0 && nrow(beta_setup) > 0) {
    beta_reorder <- match(colnames(x_matrix), beta_setup$beta_label)
    if (any(is.na(beta_reorder))) {
      warning("some beta labels not found in X matrix columns")
      beta_reorder <- beta_reorder[!is.na(beta_reorder)]
      if (length(beta_reorder) == 0) {
        beta_setup <- data.frame()
      } else {
        beta_setup <- beta_setup[beta_reorder, ]
      }
    } else {
      beta_setup <- beta_setup[beta_reorder, ]
    }
  }

  if (ncol(a_matrix) > 0 && nrow(gamma_setup) > 0) {
    gamma_reorder <- match(colnames(a_matrix), gamma_setup$gamma_label)
    if (any(is.na(gamma_reorder))) {
      warning("problem with random names")
      gamma_reorder <- gamma_reorder[!is.na(gamma_reorder)]
      if (length(gamma_reorder) == 0) {
        gamma_setup <- data.frame()
      } else {
        gamma_setup <- gamma_setup[gamma_reorder, ]
      }
    } else {
      gamma_setup <- gamma_setup[gamma_reorder, ]
    }
  }

  resp_idx <- which(vapply(terms, function(t) {
    inherits(t, "response")
  }, logical(1L)))
  obs_idx <- which(vapply(terms, function(t) {
    methods::is(t, "model") && !is.na(t@ad_kind) &&
      identical(t@ad_kind, "observations")
  }, logical(1L)))
  y <- numeric(0)
  term_idx <- integer(0)
  if (length(resp_idx) == 1L) {
    term_idx <- resp_idx[1L]
  } else if (length(resp_idx) == 0L && length(obs_idx) == 1L) {
    term_idx <- obs_idx[1L]
  } else if (length(resp_idx) > 1L || length(obs_idx) > 1L) {
    warning("multiple response or observation-density terms; using the first")
    term_idx <- if (length(resp_idx) > 0L) resp_idx[1L] else obs_idx[1L]
  }
  if (length(term_idx) == 1L) {
    y <- as.numeric(data[[terms[[term_idx]]@term]])
  } else if (methods::is(formula_in, "formula")) {
    lhs <- formula_in[[2L]]
    if (is.symbol(lhs)) {
      lhs_name <- as.character(lhs)
      if (lhs_name %in% names(data)) {
        y <- as.numeric(data[[lhs_name]])
      } else {
        warning("response variable '", lhs_name, "' not found in data")
      }
    } else {
      warning("no response variable")
    }
  } else {
    warning("no response variable")
  }

  if (is.null(beta_setup)) {
    beta_theta_names <- names(theta_setup)
  } else {
    beta_theta_names <- intersect(colnames(beta_setup), colnames(theta_setup))
  }

  beta_theta_names <- setdiff(beta_theta_names, "order")

  build_parameters_info <- function(beta_df, theta_df, names_common) {
    beta_rows <- if (!is.null(beta_df) && nrow(beta_df) > 0L) {
      out <- beta_df[, names_common, drop = FALSE]
      out$transform <- FALSE
      out
    } else {
      NULL
    }
    theta_rows <- if (nrow(theta_df) > 0L) {
      cols <- intersect(c(names_common, "transform"), names(theta_df))
      theta_df[, cols, drop = FALSE]
    } else {
      NULL
    }
    out <- do.call(rbind, Filter(Negate(is.null), list(beta_rows, theta_rows)))
    if (is.null(out) || nrow(out) == 0L) {
      return(data.frame(
        term = character(),
        model = character(),
        label = character(),
        init = numeric(),
        lower = numeric(),
        upper = numeric(),
        parscale = numeric(),
        type = factor(levels = .type_factor_levels),
        transform = logical(0),
        stringsAsFactors = FALSE
      ))
    }
    rownames(out) <- NULL
    transform <- out$transform
    transform[is.na(transform)] <- FALSE
    idx <- which(transform)
    if (length(idx) > 0L) {
      for (col in c("init", "lower", "upper")) {
        if (col %in% names(out)) {
          out[idx, col] <- log(out[idx, col])
        }
      }
    }
    out
  }

  list(
    y = y,
    A = a_matrix,
    X = x_matrix,
    data = data,
    info = list(
      beta = beta_setup,
      gamma = gamma_setup,
      theta = theta_setup,
      parameters = build_parameters_info(
        beta_setup,
        theta_setup,
        beta_theta_names
      )
    ),
    terms = terms
  )
}
