#' Set up model for adlaplace
#'
#' Function to set up model for adlaplace by extracting outcome variable and covariates,
#' creating design matrices, and preparing data structures for adlaplace optimization.
#'
#' @param formula Model formula containing model terms
#' @param data Data frame containing the variables
#' @param verbose Logical indicating whether to print verbose output (default: FALSE)
#'
#' @return An object of class \code{ad_model} (S4) from \code{\link{ad_model}}.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' model_stuff <- model_setup(
#'   data = dat,
#'   formula = y ~ x1 + f(x2, model = "iwp", p = 2, knots = seq(0, 1, len = 11)) + f(fac, model = "iid")
#' )
#' }
#'
#' @export
model_setup <- function(formula, data, verbose = FALSE) {
  # Extract outcome variable and covariates

  # Get model terms and collect information
  if (methods::is(formula, "formula")) {
    the_terms <- collect_terms(
      formula,
      verbose = verbose
    )
  } else {
    the_terms <- formula
    if (inherits(the_terms, "model")) {
      the_terms <- list(the_terms)
    }
    inherits_seq <- unlist(lapply(the_terms, inherits, what = "model"))
    if (!all(inherits_seq)) {
      warning("formula must be of class formula or a list of objects which inherit class model")
    }
  }

  the_terms = lapply(the_terms, add_by_levels, data)

  # Create design matrices
  design_list <- lapply(
    the_terms,
    design,
    data = data
  )

  # Identify fixed and random effects
  terms_with_gamma <- sapply(the_terms, methods::slot, "type") == "random"
  terms_with_beta <- sapply(the_terms, methods::slot, "type") == "fixed"

  # Build X (fixed effects) and A (random effects) matrices
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

  # Get parameter information
  theta_info_list <- Filter(Negate(is.null), lapply(the_terms, theta_info))
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

  beta_setup <- do.call(rbind, lapply(the_terms, beta_info, data = data))

  random_info_list <- Filter(Negate(is.null), lapply(the_terms, random_info, data = data))
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

  # Reorder beta and gamma setups to match matrix column names
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
      print(setdiff(gamma_setup$gamma_label, colnames(a_matrix)))
      print(setdiff(colnames(a_matrix), gamma_setup$gamma_label))
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

  # Create gamma-theta mapping matrix
  gamma_setup_sub <- gamma_setup[!is.na(gamma_setup$theta_id), ]
  if (nrow(gamma_setup_sub) > 0) {
    gamma_theta_map <- Matrix::sparseMatrix(
      i = gamma_setup_sub$id, j = gamma_setup_sub$theta_id,
      x = 1.0,
      dims = c(nrow(gamma_setup), nrow(theta_setup)),
      index1 = FALSE
    )
  } else {
    gamma_theta_map <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(0, nrow(theta_setup))
    )
  }

  gamma_map <- as_ngC(gamma_theta_map)

  # Validate column names
  if (ncol(a_matrix) > 0 && nrow(gamma_setup) > 0) {
    if (!(all(colnames(a_matrix) == gamma_setup$gamma_label))) {
      warning("some names of A don't match up")
      print(table(colnames(a_matrix) %in% gamma_setup$gamma_label))
      print(utils::str(setdiff(colnames(a_matrix), gamma_setup$gamma_label)))
      print(utils::str(setdiff(gamma_setup$gamma_label, colnames(a_matrix))))
    }
  }

  if (ncol(x_matrix) > 0 && nrow(beta_setup) > 0) {
    if (!(all(colnames(x_matrix) == beta_setup$beta_name))) {
      warning("some names of X don't match beta_setup")
      print(table(colnames(x_matrix) %in% beta_setup$beta_name))
      print(utils::str(setdiff(colnames(x_matrix), beta_setup$beta_name)))
      print(utils::str(setdiff(beta_setup$beta_name, colnames(x_matrix))))
    }
  }


  # Convert to appropriate sparse matrix types
  if (ncol(a_matrix) > 0) {
    ATp <- Matrix::t(a_matrix)
    if (!inherits(ATp, "CsparseMatrix")) {
      ATp <- methods::as(ATp, "CsparseMatrix")
    }
  } else {
    ATp <- methods::as(
      Matrix::sparseMatrix(
        i = c(), j = c(),
        dims = c(0, nrow(data))
      ),
      "dMatrix"
    )
  }

  if (ncol(x_matrix) > 0) {
    XTp <- Matrix::t(x_matrix)
    if (!inherits(XTp, "CsparseMatrix")) {
      XTp <- methods::as(XTp, "CsparseMatrix")
    }
  } else {
    XTp <- methods::as(
      Matrix::sparseMatrix(
        i = c(), j = c(),
        dims = c(0, nrow(data))
      ), "dMatrix"
    )
  }

  elgmMatrix <- as_ngC(
    Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      dims = c(nrow(data), 0L),
      giveCsparse = TRUE
    )
  )

  n_beta <- nrow(XTp)
  n_theta <- nrow(theta_setup)
  if (n_beta > 0L) {
    beta_map <- as_ngC(Matrix::Diagonal(n_beta))
  } else {
    beta_map <- empty_ngC()
  }
  if (n_theta > 0L) {
    theta_map <- as_ngC(Matrix::Diagonal(n_theta))
  } else {
    theta_map <- empty_ngC()
  }

  the_response <- which(
    unlist(lapply(the_terms, function(xx) any(class(xx) == "response")))
  )
  y <- numeric(0)
  if (length(the_response) == 1) {
    y <- as.numeric(data[[
      the_terms[[
        the_response[1]
      ]]@term
    ]])
  } else {
    warning("no response variable")
  }

  if (is.null(beta_setup)) {
    beta_theta_names <- names(theta_setup)
  } else {
    beta_theta_names <- intersect(colnames(beta_setup), colnames(theta_setup))
  }

  beta_theta_names <- setdiff(beta_theta_names, "order")

  ad_model(
    y = y,
    ATp = ATp,
    XTp = XTp,
    beta_map = beta_map,
    gamma_map = gamma_map,
    theta_map = theta_map,
    elgmMatrix = elgmMatrix,
    data = data,
    terms = the_terms,
    info = list(
      beta = beta_setup,
      gamma = gamma_setup,
      theta = theta_setup,
      parameters = rbind(
        beta_setup[, beta_theta_names, drop = FALSE],
        theta_setup[, beta_theta_names, drop = FALSE]
      )
    )
  )
}
