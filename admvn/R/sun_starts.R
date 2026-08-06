#' Map mean and covariance to a SUN(2,2) start (zero / tiny skewness)
#'
#' Uses marginal standard deviations and unconstrained correlation free
#' parameters for \code{(nu, omega)}. Skew loadings and residual-correlation
#' free parameters are zero, or a tiny diagonal loading if
#' \code{skew_strength > 0}.
#'
#' @param mu Length-2 mean.
#' @param Sigma \eqn{2\times 2} covariance.
#' @param skew_strength Diagonal loading for \code{L*} (default 0).
#' @return Named length-10 parameter vector for [make_sun22_params()].
#' @export
make_sun22_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 2L) stop("mu must have length 2")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(2L, 2L))) stop("Sigma must be 2 x 2")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  c(
    xi1 = mu[1], xi2 = mu[2],
    nu1 = nu[1], nu2 = nu[2],
    omega12 = om[1],
    L11 = skew_strength, L12 = 0,
    L21 = 0, L22 = skew_strength,
    gamma12 = 0
  )
}

#' SUN(2,2) start from weighted sample moments on a quadrature
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun22_start_from_normal()].
#' @return Named length-10 parameter vector.
#' @export
make_sun22_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun22_start_from_normal(mu, S, skew_strength = skew_strength)
}

#' Map mean and covariance to a SUN(3,3) start (zero / tiny skewness)
#'
#' Uses marginal standard deviations and unconstrained correlation free
#' parameters for \code{(nu, omega)}. Skew loadings and residual-correlation
#' free parameters are zero, or a tiny diagonal loading if
#' \code{skew_strength > 0}.
#'
#' @param mu Length-3 mean.
#' @param Sigma \eqn{3\times 3} covariance.
#' @param skew_strength Diagonal loading for \code{L*} (default 0).
#' @return Named length-21 parameter vector for [make_sun_params()].
#' @export
make_sun33_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 3L) stop("mu must have length 3")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(3L, 3L))) stop("Sigma must be 3 x 3")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  c(
    xi1 = mu[1], xi2 = mu[2], xi3 = mu[3],
    nu1 = nu[1], nu2 = nu[2], nu3 = nu[3],
    omega12 = om[1], omega13 = om[2], omega23 = om[3],
    L11 = skew_strength, L12 = 0, L13 = 0,
    L21 = 0, L22 = skew_strength, L23 = 0,
    L31 = 0, L32 = 0, L33 = skew_strength,
    gamma12 = 0, gamma13 = 0, gamma23 = 0
  )
}

#' SUN(3,3) start from weighted sample moments on a quadrature
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun33_start_from_normal()].
#' @return Named length-21 parameter vector.
#' @export
make_sun33_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun33_start_from_normal(mu, S, skew_strength = skew_strength)
}

#' Map mean and covariance to a SUN(4,4) start (zero / tiny skewness)
#'
#' @param mu Length-4 mean.
#' @param Sigma \eqn{4\times 4} covariance.
#' @param skew_strength Diagonal loading for \code{L*} (default 0).
#' @return Named length-36 parameter vector for [make_sun44_params()].
#' @export
make_sun44_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 4L) stop("mu must have length 4")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(4L, 4L))) stop("Sigma must be 4 x 4")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  c(
    xi1 = mu[1], xi2 = mu[2], xi3 = mu[3], xi4 = mu[4],
    nu1 = nu[1], nu2 = nu[2], nu3 = nu[3], nu4 = nu[4],
    omega12 = om[1], omega13 = om[2], omega23 = om[3],
    omega14 = om[4], omega24 = om[5], omega34 = om[6],
    L11 = skew_strength, L12 = 0, L13 = 0, L14 = 0,
    L21 = 0, L22 = skew_strength, L23 = 0, L24 = 0,
    L31 = 0, L32 = 0, L33 = skew_strength, L34 = 0,
    L41 = 0, L42 = 0, L43 = 0, L44 = skew_strength,
    gamma12 = 0, gamma13 = 0, gamma23 = 0,
    gamma14 = 0, gamma24 = 0, gamma34 = 0
  )
}

#' SUN(4,4) start from weighted sample moments on a quadrature
#'
#' @inheritParams make_sun33_start_from_gh
#' @return Named length-36 parameter vector.
#' @export
make_sun44_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun44_start_from_normal(mu, S, skew_strength = skew_strength)
}
