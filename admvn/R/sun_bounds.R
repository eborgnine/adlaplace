#' Default L-BFGS-B box constraints for SUN(2,2) parameters
#'
#' Diagonal loadings \code{L11,L22} are in \code{[-1, 1]} (correlation
#' scale). Off-diagonal loadings \code{L12,L21} are in \code{[-0.5, 0.5]},
#' matching [sun33_bounds()].
#'
#' @return List with \code{lower} and \code{upper} length-10 numeric vectors.
#' @seealso [make_sun22_params()], [fit_sun22_quad()]
#' @export
sun22_bounds <- function() {
  list(
    lower = c(
      -Inf, -Inf,
      1e-6, 1e-6,
      -Inf,
      -1, -0.5, -0.5, -1, # L11, L12, L21, L22
      -5
    ),
    upper = c(
      Inf, Inf,
      Inf, Inf,
      Inf,
      1, 0.5, 0.5, 1,
      5
    )
  )
}

#' Default box constraints for SUN(2,2) hyperspherical parameters
#'
#' Independent-SN pair slots use \code{atanh(+/-0.8)} and all other joint
#' coordinates use \code{atanh(+/-0.5)}.
#' @return List of length-10 \code{lower} and \code{upper} vectors.
#' @export
sun22_hs_bounds <- function() {
  z_lo <- rep(atanh(-0.5), 6L)
  z_hi <- rep(atanh(0.5), 6L)
  pair <- .sun_hs_pair_z_idx(2L)
  z_lo[pair] <- atanh(-0.8)
  z_hi[pair] <- atanh(0.8)
  list(
    lower = c(rep(-Inf, 2L), rep(1e-6, 2L), z_lo),
    upper = c(rep(Inf, 4L), z_hi)
  )
}

#' Default L-BFGS-B box constraints for SUN(3,3) parameters
#'
#' Diagonal loadings \code{L11,L22,L33} are in \code{[-1, 1]} (correlation
#' scale). Off-diagonal loadings \code{Lij}/\code{Lji} are in
#' \code{[-0.5, 0.5]}.
#'
#' @return List with \code{lower} and \code{upper} length-21 numeric vectors.
#' @seealso [make_sun_params()], [fit_sun33_quad()]
#' @export
sun33_bounds <- function() {
  list(
    lower = c(
      -Inf, -Inf, -Inf,
      1e-6, 1e-6, 1e-6,
      -Inf, -Inf, -Inf,
      -1, -0.5, -0.5, # L11, L12, L13
      -0.5, -1, -0.5, # L21, L22, L23
      -0.5, -0.5, -1, # L31, L32, L33
      -5, -5, -5
    ),
    upper = c(
      Inf, Inf, Inf,
      Inf, Inf, Inf,
      Inf, Inf, Inf,
      1, 0.5, 0.5,
      0.5, 1, 0.5,
      0.5, 0.5, 1,
      5, 5, 5
    )
  )
}

#' Default L-BFGS-B box constraints for SUN(3,3) hyperspherical parameters
#'
#' Soft boxes on the joint unit-row Cholesky coordinates of [make_sun_hs_params()]:
#' independent-SN pair slots \code{z41,z52,z63} lie in
#' \code{atanh(+/-0.8)}; all other \code{z*} lie in \code{atanh(+/-0.5)}.
#' These are regularization / start boxes, not a Gershgorin PD certificate.
#'
#' @return List with \code{lower} and \code{upper} length-21 numeric vectors.
#' @seealso [make_sun_hs_params()], [dsun_hs()]
#' @export
sun33_hs_bounds <- function() {
  z_lo <- rep(atanh(-0.5), 15L)
  z_hi <- rep(atanh(0.5), 15L)
  pair <- .sun_hs_pair_z_idx(3L)
  z_lo[pair] <- atanh(-0.8)
  z_hi[pair] <- atanh(0.8)
  list(
    lower = c(
      -Inf, -Inf, -Inf,
      1e-6, 1e-6, 1e-6,
      z_lo
    ),
    upper = c(
      Inf, Inf, Inf,
      Inf, Inf, Inf,
      z_hi
    )
  )
}

#' Default L-BFGS-B box constraints for SUN(4,4) parameters
#'
#' Skew loadings \code{Lij} act on the correlation scale, so they are
#' restricted to \code{[-1, 1]}.
#'
#' @return List with \code{lower} and \code{upper} length-36 numeric vectors.
#' @seealso [make_sun44_params()], [fit_sun44_quad()]
#' @export
sun44_bounds <- function() {
  list(
    lower = c(
      rep(-Inf, 4),
      rep(1e-6, 4),
      rep(-Inf, 6),
      rep(-1, 16),
      rep(-5, 6)
    ),
    upper = c(
      rep(Inf, 4),
      rep(Inf, 4),
      rep(Inf, 6),
      rep(1, 16),
      rep(5, 6)
    )
  )
}

#' Default box constraints for SUN(4,4) hyperspherical parameters
#'
#' Independent-SN pair slots use \code{atanh(+/-0.8)} and all other joint
#' coordinates use \code{atanh(+/-0.5)}.
#' @return List of length-36 \code{lower} and \code{upper} vectors.
#' @export
sun44_hs_bounds <- function() {
  z_lo <- rep(atanh(-0.5), 28L)
  z_hi <- rep(atanh(0.5), 28L)
  pair <- .sun_hs_pair_z_idx(4L)
  z_lo[pair] <- atanh(-0.8)
  z_hi[pair] <- atanh(0.8)
  list(
    lower = c(rep(-Inf, 4L), rep(1e-6, 4L), z_lo),
    upper = c(rep(Inf, 8L), z_hi)
  )
}
