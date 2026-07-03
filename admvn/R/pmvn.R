#' Multivariate normal CDF with derivatives
#'
#' Evaluates \eqn{P(\code{lower} < X < \code{upper})} for
#' \eqn{X \sim N(\code{mean}, \code{sigma})} using Genz's quasi-Monte Carlo
#' algorithm. Returns the probability, a Monte Carlo error estimate, and
#' derivatives of the approximation with respect to \code{upper}.
#'
#' @param upper Numeric vector of upper integration limits.
#' @param lower Numeric vector of lower integration limits. Defaults to
#'   \code{-Inf}.
#' @param mean Numeric mean vector. Defaults to zero.
#' @param sigma Positive-definite covariance matrix.
#' @param n_points Number of quasi-Monte Carlo points per random shift.
#' @param n_shifts Number of random shifts for the QMC rule.
#' @param seed Random seed for generating QMC shifts.
#'
#' @return A list with components \code{value}, \code{error}, \code{gradient},
#'   and \code{hessian}. The gradient and Hessian are with respect to
#'   \code{upper}.
#'
#' @examples
#' sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' pmvn(c(1, 1), lower = c(-Inf, -Inf), sigma = sigma)
#'
#' @export
pmvn <- function(upper,
                 lower = -Inf,
                 mean = 0,
                 sigma,
                 n_points = 1021L,
                 n_shifts = 8L,
                 seed = 1L) {
  upper <- as.numeric(upper)
  n <- length(upper)
  if (length(lower) == 1L) {
    lower <- rep(lower, n)
  } else {
    lower <- as.numeric(lower)
  }
  if (length(mean) == 1L) {
    mean <- rep(mean, n)
  } else {
    mean <- as.numeric(mean)
  }
  sigma <- as.matrix(sigma)
  if (nrow(sigma) != n || ncol(sigma) != n) {
    stop("sigma must be ", n, " x ", n)
  }
  if (length(lower) != n || length(mean) != n) {
    stop("lower and mean must have length ", n, " or 1")
  }
  if (any(lower > upper)) {
    stop("lower must be <= upper for all components")
  }
  if (n == 1L) {
    sd <- sqrt(sigma[1, 1])
    z_u <- (upper - mean) / sd
    z_l <- if (is.infinite(lower) && lower < 0) -Inf else (lower - mean) / sd
    value <- pnorm(z_u) - pnorm(z_l)
    grad <- dnorm(z_u) / sd
    hess <- matrix(-z_u * grad / sd, 1L, 1L)
    return(list(
      value = value,
      error = 0,
      gradient = grad,
      hessian = hess
    ))
  }
  storage.mode(sigma) <- "double"
  chol(sigma)
  pmvn_cpp(
    upper = upper,
    lower = lower,
    mean = mean,
    sigma = sigma,
    n_points = as.integer(n_points),
    n_shifts = as.integer(n_shifts),
    seed = as.integer(seed)
  )
}
