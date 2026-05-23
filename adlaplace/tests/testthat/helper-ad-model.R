#' Build a minimal ad_model for shard tests
#' @keywords internal
test_ad_model <- function(y, ATp, XTp, config, theta_local_row = 0L) {
  adlaplace:::ad_model_from_config_matrices(
    y = y,
    ATp = ATp,
    XTp = XTp,
    config = config,
    theta_local_row = theta_local_row
  )
}

#' Random term model + precision for tests using a single Amat block
#' @keywords internal
test_random_shard <- function(model, config, gamma_ids, theta_id, Q) {
  n_gamma_full <- length(config$gamma)
  gamma_map <- adlaplace:::as_ngC(
    Matrix::sparseMatrix(
      i = as.integer(gamma_ids),
      j = rep(0L, length(gamma_ids)),
      x = rep(1, length(gamma_ids)),
      dims = c(n_gamma_full, 1L),
      index1 = FALSE,
      giveCsparse = TRUE
    )
  )
  term_model <- adlaplace:::ad_model_random_maps(
    gamma_map = gamma_map,
    theta_map = model@theta_map,
    n_beta = length(config$beta),
    n_theta = length(config$theta)
  )
  list(
    model = term_model,
    precision = list(Q = Q)
  )
}
