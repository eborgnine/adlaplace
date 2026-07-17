#' Build a minimal ad_data for shard tests
#'
#' Optional \code{ad_kind} / \code{ad_fun} / \code{precision} are stored on
#' the returned object so it can be passed directly to \code{ad_fun_ptr()}.
#' @keywords internal
test_ad_data <- function(y, A, X, config,
                         theta_local_row = 0L,
                         ad_kind = NA_character_,
                         ad_fun = NA_character_,
                         precision = NULL) {
  data <- adlaplace:::ad_data_from_config_matrices(
    y = y,
    A = A,
    X = X,
    config = config,
    theta_local_row = theta_local_row
  )
  data@ad_kind <- ad_kind
  data@ad_fun <- ad_fun
  data@precision <- precision
  data
}

#' Set \code{ad_kind}, \code{ad_fun} and \code{precision} on an \code{ad_data}
#' so it can be passed to \code{ad_fun_ptr()}.
#' @keywords internal
as_shard <- function(data, ad_kind, ad_fun, precision = NULL) {
  data@ad_kind <- ad_kind
  data@ad_fun <- ad_fun
  data@precision <- precision
  data
}

#' Random-shard ad_data for tests with a single Amat block
#' @keywords internal
test_random_shard <- function(data, config, gamma_ids, theta_id, Q,
                              ad_fun = "random_diagonal") {
  n_gamma_full <- length(config$gamma)
  n_term <- length(gamma_ids)
  gamma_map <- Matrix::sparseMatrix(
    i = as.integer(gamma_ids),
    j = seq.int(0L, length.out = n_term),
    dims = c(n_gamma_full, n_term),
    index1 = FALSE,
    repr = "C"
  )
  adlaplace:::ad_data(
    beta_map = Matrix::Matrix(nrow = length(config$beta), ncol = 0L),
    gamma_map = gamma_map,
    theta_map = data@theta_map,
    ad_kind = "random",
    ad_fun = ad_fun,
    precision = Q
  )
}
