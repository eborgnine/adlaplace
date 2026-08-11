#' Build SUN(3,2) parameters from a joint hyperspherical vector
#'
#' The joint correlation of \eqn{(U,V)\in\mathbb{R}^5} is formed from ten
#' unconstrained unit-row Cholesky coordinates in row-major order.
#'
#' @param par Numeric vector of length 16: \code{xi(3)}, \code{nu(3)}, then
#'   \code{z21,z31,z32,z41,z42,z43,z51,z52,z53,z54}.
#' @return A SUN \code{dp} list with \code{xi}, \code{Omega} (\eqn{3\times 3}),
#'   \code{Delta} (\eqn{3\times 2}), zero \code{tau} (length 2), and
#'   \code{Gamma} (\eqn{2\times 2}).
#' @seealso [dsun32_hs()], [sun32_hs_bounds()]
#' @export
make_sun32_hs_params <- function(par) {
  if (length(par) != 16L) {
    stop(
      "par must have length 16: xi(3), nu(3), z(10) in order ",
      paste(.sun_hs_z_names(3L, 2L), collapse = ", "),
      call. = FALSE
    )
  }
  z_names <- .sun_hs_z_names(3L, 2L)
  names(par) <- c("xi1", "xi2", "xi3", "nu1", "nu2", "nu3", z_names)
  xi <- unname(par[c("xi1", "xi2", "xi3")])
  nu <- unname(par[c("nu1", "nu2", "nu3")])
  J <- .corr_from_free(unname(par[z_names]), 5L, eps = 1e-6)
  Cu <- J[1:3, 1:3, drop = FALSE]
  list(
    xi = xi,
    Omega = unname(outer(pmax(nu, 1e-12), pmax(nu, 1e-12)) * Cu),
    Delta = unname(J[1:3, 4:5, drop = FALSE]),
    tau = rep(0, 2L),
    Gamma = unname(J[4:5, 4:5, drop = FALSE])
  )
}
