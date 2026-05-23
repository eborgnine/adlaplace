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
  parent_layout <- sizes(model)

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
  gamma_map <- Matrix::sparseMatrix(
    i = parent_gamma_id,
    j = seq.int(0L, length.out = n_term),
    dims = c(parent_layout$n_gamma, n_term),
    index1 = FALSE,
    giveCsparse = TRUE
  )

  ad_model(
    beta_map = Matrix::Matrix(nrow = parent_layout$n_beta, ncol = 0L),
    gamma_map = gamma_map,
    theta_map = model@theta_map
  )
}
