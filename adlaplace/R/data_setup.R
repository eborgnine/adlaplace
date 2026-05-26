#' Build design matrices and parameter metadata from a formula
#'
#' @param formula Model formula or named list of term objects.
#' @param data Data frame.
#' @param verbose Print parsing messages.
#' @return List with \code{terms}, \code{y}, \code{A}, \code{X}, \code{data}, \code{info}.
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
      stringsAsFactors = FALSE
    )
  } else {
    theta_setup <- do.call(rbind, theta_info_list)
    theta_setup <- theta_setup[order(theta_setup$type, theta_setup$label), ]
  }

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

  list(
    y = y,
    A = a_matrix,
    X = x_matrix,
    data = data,
    info = list(
      beta = beta_setup,
      gamma = gamma_setup,
      theta = theta_setup,
      parameters = rbind(
        beta_setup[, beta_theta_names, drop = FALSE],
        theta_setup[, beta_theta_names, drop = FALSE]
      )
    ),
    terms = terms
  )
}
