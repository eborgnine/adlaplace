#' Gauss-Hermite nodes and weights under a multivariate normal
#'
#' Product Gauss-Hermite quadrature for \eqn{N(\mu,\Sigma)}, returned as
#' abscissae in the original scale and corresponding product weights.
#'
#' @param mu Mean vector.
#' @param Sigma Covariance matrix (made symmetric; Cholesky used for the map).
#' @param n Number of univariate GH nodes per dimension (default 15).
#'
#' @return A list with \code{x} (\eqn{n^d \times d} matrix) and \code{w}
#'   (length \eqn{n^d} weights).
#'
#' @export
gh_points <- function(mu, Sigma, n = 15) {
  mu <- as.numeric(mu)
  d <- length(mu)
  gh <- statmod::gauss.quad.prob(n, dist = "normal")
  ids <- expand.grid(rep(list(seq_len(n)), d))

  X0 <- as.matrix(
    data.frame(lapply(seq_len(d), function(k) gh$nodes[ids[, k]]))
  )
  w0 <- apply(ids, 1, function(ii) prod(gh$weights[ii]))

  Sigma <- (as.matrix(Sigma) + t(as.matrix(Sigma))) / 2
  X <- sweep(X0 %*% chol(Sigma), 2, mu, "+")
  list(x = X, w = w0)
}

#' Weighted column means
#'
#' @param X Numeric matrix of observations (rows).
#' @param w Nonnegative weights (recycled / renormalized).
#' @return Named numeric vector of length \code{ncol(X)}.
#' @export
weighted_mean <- function(X, w) {
  X <- as.matrix(X)
  w <- as.numeric(w)
  w <- w / sum(w)
  colSums(X * w)
}

#' Weighted covariance
#'
#' @inheritParams weighted_mean
#' @return Symmetric \code{ncol(X)} covariance matrix.
#' @export
weighted_cov <- function(X, w) {
  X <- as.matrix(X)
  w <- as.numeric(w)
  w <- w / sum(w)
  mu <- weighted_mean(X, w)
  Xc <- sweep(X, 2, mu)
  crossprod(Xc * sqrt(w), Xc * sqrt(w))
}
