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

#' Map mean and covariance to a SUN(2,2) hyperspherical start
#' @inheritParams make_sun22_start_from_normal
#' @return Named length-10 vector for [make_sun22_hs_params()].
#' @export
make_sun22_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 2L) stop("mu must have length 2")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(2L, 2L))) stop("Sigma must be 2 x 2")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  z <- setNames(numeric(6L), .sun_hs_z_names(2L, 2L))
  z["z21"] <- .corr_free_from_C(stats::cov2cor(Sigma))
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z31", "z42")] <- atanh(s)
  }
  c(xi1 = mu[1], xi2 = mu[2], nu1 = nu[1], nu2 = nu[2], z)
}

#' SUN(2,2) hyperspherical start from weighted sample moments
#' @inheritParams make_sun22_start_from_gh
#' @return Named length-10 parameter vector.
#' @export
make_sun22_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  make_sun22_hs_start_from_normal(
    weighted_mean(Z, w), weighted_cov(Z, w), skew_strength)
}

#' Map mean and covariance to a SUN(3,2) hyperspherical start
#'
#' Exact Gaussian start on the joint hyperspherical layout of
#' [make_sun32_hs_params()]: the three \eqn{U}-block free coordinates come from
#' inverting \code{cov2cor(Sigma)}; all cross / \eqn{V} coordinates are zero.
#' Optional \code{skew_strength} sets pair slots \code{z41,z52} to
#' \code{atanh(skew_strength)}.
#'
#' @param mu Length-3 mean.
#' @param Sigma \eqn{3\times 3} covariance.
#' @param skew_strength Pair-slot \eqn{t=\tanh(z)} value in \eqn{(-1,1)}
#'   (default 0).
#' @return Named length-16 parameter vector for [make_sun32_hs_params()].
#' @export
make_sun32_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 3L) stop("mu must have length 3")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(3L, 3L))) stop("Sigma must be 3 x 3")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  z <- setNames(numeric(10L), .sun_hs_z_names(3L, 2L))
  z[c("z21", "z31", "z32")] <- om
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z41", "z52")] <- atanh(s)
  }

  c(
    xi1 = mu[1], xi2 = mu[2], xi3 = mu[3],
    nu1 = nu[1], nu2 = nu[2], nu3 = nu[3],
    z
  )
}

#' SUN(3,2) hyperspherical start from weighted sample moments
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun32_hs_start_from_normal()].
#' @return Named length-16 parameter vector.
#' @export
make_sun32_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun32_hs_start_from_normal(mu, S, skew_strength = skew_strength)
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

#' Map mean and covariance to a SUN(3,3) hyperspherical start
#'
#' Exact Gaussian start on the joint hyperspherical layout of
#' [make_sun_hs_params()]: the three \eqn{U}-block free coordinates come from
#' inverting \code{cov2cor(Sigma)}; all cross / \eqn{V} coordinates are zero
#' (so \eqn{J=\operatorname{blockdiag}(C_u,I)}). Optional
#' \code{skew_strength} sets the independent-SN pair slots to
#' \code{atanh(skew_strength)}.
#'
#' @param mu Length-3 mean.
#' @param Sigma \eqn{3\times 3} covariance.
#' @param skew_strength Pair-slot \eqn{t=\tanh(z)} value in \eqn{(-1,1)}
#'   (default 0).
#' @return Named length-21 parameter vector for [make_sun_hs_params()].
#' @export
make_sun33_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 3L) stop("mu must have length 3")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(3L, 3L))) stop("Sigma must be 3 x 3")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  z <- setNames(numeric(15L), .sun_hs_z_names(3L, 3L))
  z[c("z21", "z31", "z32")] <- om
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z41", "z52", "z63")] <- atanh(s)
  }

  c(
    xi1 = mu[1], xi2 = mu[2], xi3 = mu[3],
    nu1 = nu[1], nu2 = nu[2], nu3 = nu[3],
    z
  )
}

#' SUN(3,3) hyperspherical start from weighted sample moments
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun33_hs_start_from_normal()].
#' @return Named length-21 parameter vector.
#' @export
make_sun33_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun33_hs_start_from_normal(mu, S, skew_strength = skew_strength)
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

#' Map mean and covariance to a SUN(4,4) hyperspherical start
#' @inheritParams make_sun44_start_from_normal
#' @return Named length-36 vector for [make_sun44_hs_params()].
#' @export
make_sun44_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 4L) stop("mu must have length 4")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(4L, 4L))) stop("Sigma must be 4 x 4")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  z <- setNames(numeric(28L), .sun_hs_z_names(4L, 4L))
  z[seq_len(6L)] <- .corr_free_from_C(stats::cov2cor(Sigma))
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z51", "z62", "z73", "z84")] <- atanh(s)
  }
  c(setNames(mu, paste0("xi", 1:4)), setNames(nu, paste0("nu", 1:4)), z)
}

#' SUN(4,4) hyperspherical start from weighted sample moments
#' @inheritParams make_sun44_start_from_gh
#' @return Named length-36 parameter vector.
#' @export
make_sun44_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  make_sun44_hs_start_from_normal(
    weighted_mean(Z, w), weighted_cov(Z, w), skew_strength)
}


#' Map mean and covariance to a SUN(4,3) hyperspherical start
#'
#' Exact Gaussian start on the joint hyperspherical layout of
#' [make_sun43_hs_params()]: the six \eqn{U}-block free coordinates come from
#' inverting \code{cov2cor(Sigma)}; all cross / \eqn{V} coordinates are zero.
#' Optional \code{skew_strength} sets pair slots \code{z51,z62,z73} to
#' \code{atanh(skew_strength)}.
#'
#' @param mu Length-4 mean.
#' @param Sigma \eqn{4\times 4} covariance.
#' @param skew_strength Pair-slot \eqn{t=\tanh(z)} value in \eqn{(-1,1)}
#'   (default 0).
#' @return Named length-29 parameter vector for [make_sun43_hs_params()].
#' @export
make_sun43_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 4L) stop("mu must have length 4")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(4L, 4L))) stop("Sigma must be 4 x 4")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  z <- setNames(numeric(21L), .sun_hs_z_names(4L, 3L))
  z[c("z21", "z31", "z32", "z41", "z42", "z43")] <- om
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z51", "z62", "z73")] <- atanh(s)
  }

  c(
    setNames(mu, paste0("xi", 1:4)),
    setNames(nu, paste0("nu", 1:4)),
    z
  )
}

#' SUN(4,3) hyperspherical start from weighted sample moments
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun43_hs_start_from_normal()].
#' @return Named length-29 parameter vector.
#' @export
make_sun43_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun43_hs_start_from_normal(mu, S, skew_strength = skew_strength)
}

#' Map mean and covariance to a SUN(4,2) hyperspherical start
#'
#' Exact Gaussian start on the joint hyperspherical layout of
#' [make_sun42_hs_params()]: the six \eqn{U}-block free coordinates come from
#' inverting \code{cov2cor(Sigma)}; all cross / \eqn{V} coordinates are zero.
#' Optional \code{skew_strength} sets pair slots \code{z51,z62} to
#' \code{atanh(skew_strength)}.
#'
#' @param mu Length-4 mean.
#' @param Sigma \eqn{4\times 4} covariance.
#' @param skew_strength Pair-slot \eqn{t=\tanh(z)} value in \eqn{(-1,1)}
#'   (default 0).
#' @return Named length-23 parameter vector for [make_sun42_hs_params()].
#' @export
make_sun42_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 4L) stop("mu must have length 4")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(4L, 4L))) stop("Sigma must be 4 x 4")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  z <- setNames(numeric(15L), .sun_hs_z_names(4L, 2L))
  z[c("z21", "z31", "z32", "z41", "z42", "z43")] <- om
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z51", "z62")] <- atanh(s)
  }

  c(
    setNames(mu, paste0("xi", 1:4)),
    setNames(nu, paste0("nu", 1:4)),
    z
  )
}

#' SUN(4,2) hyperspherical start from weighted sample moments
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun42_hs_start_from_normal()].
#' @return Named length-23 parameter vector.
#' @export
make_sun42_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun42_hs_start_from_normal(mu, S, skew_strength = skew_strength)
}

#' Map mean and covariance to a SUN(5,2) hyperspherical start
#'
#' Exact Gaussian start on the joint hyperspherical layout of
#' [make_sun52_hs_params()]: the ten \eqn{U}-block free coordinates come from
#' inverting \code{cov2cor(Sigma)}; all cross / \eqn{V} coordinates are zero.
#' Optional \code{skew_strength} sets pair slots \code{z61,z72} to
#' \code{atanh(skew_strength)}.
#'
#' @param mu Length-5 mean.
#' @param Sigma \eqn{5\times 5} covariance.
#' @param skew_strength Pair-slot \eqn{t=\tanh(z)} value in \eqn{(-1,1)}
#'   (default 0).
#' @return Named length-31 parameter vector for [make_sun52_hs_params()].
#' @export
make_sun52_hs_start_from_normal <- function(mu, Sigma, skew_strength = 0) {
  mu <- as.numeric(mu)
  if (length(mu) != 5L) stop("mu must have length 5")
  Sigma <- as.matrix(Sigma)
  if (!all(dim(Sigma) == c(5L, 5L))) stop("Sigma must be 5 x 5")
  Sigma <- (Sigma + t(Sigma)) / 2
  nu <- sqrt(pmax(diag(Sigma), 1e-8))
  C <- stats::cov2cor(Sigma)
  om <- .corr_free_from_C(C)

  z <- setNames(numeric(21L), .sun_hs_z_names(5L, 2L))
  z[seq_len(10L)] <- om
  if (abs(skew_strength) > 0) {
    s <- max(min(as.numeric(skew_strength), 1 - 1e-8), -1 + 1e-8)
    z[c("z61", "z72")] <- atanh(s)
  }

  c(
    setNames(mu, paste0("xi", 1:5)),
    setNames(nu, paste0("nu", 1:5)),
    z
  )
}

#' SUN(5,2) hyperspherical start from weighted sample moments
#'
#' @param quad List with \code{x} and \code{w}.
#' @param skew_strength Passed to [make_sun52_hs_start_from_normal()].
#' @return Named length-31 parameter vector.
#' @export
make_sun52_hs_start_from_gh <- function(quad, skew_strength = 0.01) {
  Z <- as.matrix(quad$x)
  w <- as.numeric(quad$w)
  w <- w / sum(w)
  mu <- weighted_mean(Z, w)
  S <- weighted_cov(Z, w)
  make_sun52_hs_start_from_normal(mu, S, skew_strength = skew_strength)
}
