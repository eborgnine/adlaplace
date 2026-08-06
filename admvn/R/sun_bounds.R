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
