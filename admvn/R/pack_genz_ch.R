#' Pack covariance into Genz chol and scales
#'
#' Converts a covariance matrix into the Genz reordered correlation Cholesky
#' factor and marginal scales used by the CppAD tape. The permutation \code{perm}
#' must match the frozen order from [pmvn_fun()].
#'
#' @param sigma Positive-definite covariance matrix.
#' @param perm Integer permutation vector from a tape object.
#'
#' @return A list with components \code{scale} and \code{ch}.
#'
#' @export
pack_genz_ch <- function(sigma, perm) {
  sigma <- as.matrix(sigma)
  storage.mode(sigma) <- "double"
  chol(sigma)
  pack_genz_ch_cpp(sigma, as.integer(perm))
}
