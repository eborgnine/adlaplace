#' Precision list for a random-effect AD shard
#'
#' @param term Model term with \code{as_kind = "random"}.
#' @param model An \code{ad_model} object from \code{model_setup()}.
#' @keywords internal
random_shard_precision <- function(term, model) {
  prec <- precision(term, data = model@data)
  if (is.null(prec) || length(prec) == 0) {
    return(NULL)
  }
  list(Q = as.numeric(Matrix::diag(prec)))
}

#' Maps-only \code{ad_model} for one random-effect term
#'
#' @param term Model term with \code{as_kind = "random"}.
#' @param model An \code{ad_model} object from \code{model_setup()}.
#' @keywords internal
random_term_model <- function(term, model) {
  gamma_setup <- model@info$gamma
  theta_setup <- model@info$theta
  parent_layout <- ad_model_layout(model)

  prec <- precision(term, data = model@data)
  if (is.null(prec) || length(prec) == 0) {
    return(NULL)
  }

  gamma_cols <- colnames(prec)
  if (is.null(gamma_cols) && !is.null(rownames(prec))) {
    gamma_cols <- rownames(prec)
  }
  if (is.null(gamma_cols)) {
    warning("precision for term ", term@term, " has no dimnames; skipping shard")
    return(NULL)
  }

  gamma_match <- match(gamma_cols, gamma_setup$gamma_label)
  if (anyNA(gamma_match)) {
    warning(
      "gamma labels for term ", term@term,
      " do not match gamma_setup; dropping unmatched columns"
    )
    keep <- !is.na(gamma_match)
    gamma_match <- gamma_match[keep]
    if (length(gamma_match) == 0) {
      return(NULL)
    }
  }

  n_term <- length(gamma_match)
  parent_gamma_id <- as.integer(gamma_setup$id[gamma_match])
  gamma_map <- methods::as(
    Matrix::sparseMatrix(
      i = parent_gamma_id,
      j = rep(0L, n_term),
      x = rep(1, n_term),
      dims = c(parent_layout$n_gamma, 1L),
      index1 = FALSE,
      giveCsparse = TRUE
    ),
    "ngCMatrix"
  )

  th <- theta_info(term)
  if (!is.null(th) && nrow(th) > 0) {
    theta_row <- match(th$label[1], theta_setup$label)
    if (is.na(theta_row)) {
      stop("theta label for term ", term@term, " not found in theta_setup")
    }
    theta_local <- as.integer(theta_setup$id[theta_row])
  } else {
    theta_ids <- unique(gamma_setup$theta_id[gamma_match])
    theta_ids <- theta_ids[!is.na(theta_ids)]
    if (length(theta_ids) == 0) {
      stop("no theta index for random term ", term@term)
    }
    theta_local <- as.integer(theta_ids[1])
  }

  theta_map <- methods::as(
    model@theta_map,
    "ngCMatrix"
  )

  ad_model_random_maps(
    gamma_map = gamma_map,
    theta_map = theta_map,
    n_beta = parent_layout$n_beta,
    n_theta = parent_layout$n_theta
  )
}
