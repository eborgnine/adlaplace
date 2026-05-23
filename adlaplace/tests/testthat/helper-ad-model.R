#' Build a minimal ad_model for shard tests
#' @keywords internal
test_ad_model <- function(y, A, X, config, theta_local_row = 0L) {
  adlaplace:::ad_model_from_config_matrices(
    y = y,
    A = A,
    X = X,
    config = config,
    theta_local_row = theta_local_row
  )
}

#' Random term model + precision for tests using a single Amat block
#' @keywords internal
test_random_shard <- function(model, config, gamma_ids, theta_id, Q) {
  n_gamma_full <- length(config$gamma)
  n_term <- length(gamma_ids)
  gamma_map <- Matrix::sparseMatrix(
    i = as.integer(gamma_ids),
    j = seq.int(0L, length.out = n_term),
    dims = c(n_gamma_full, n_term),
    index1 = FALSE,
    giveCsparse = TRUE
  )
  term_model <- adlaplace::ad_model(
    beta_map = Matrix::Matrix(nrow = length(config$beta), ncol = 0L),
    gamma_map = gamma_map,
    theta_map = model@theta_map
  )
  list(
    model = term_model,
    precision = list(Q = Q)
  )
}
