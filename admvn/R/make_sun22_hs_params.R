#' Build SUN(2,2) parameters from a joint hyperspherical vector
#'
#' The joint correlation of \eqn{(U,V)} is formed from six unconstrained
#' unit-row Cholesky coordinates in row-major order.
#'
#' @param par Numeric vector of length 10: \code{xi(2)}, \code{nu(2)}, then
#'   \code{z21,z31,z32,z41,z42,z43}.
#' @return A SUN \code{dp} list with \code{xi}, \code{Omega}, \code{Delta},
#'   zero \code{tau}, and \code{Gamma}.
#' @seealso [dsun22_hs()], [make_sun22_params()], [sun22_hs_bounds()]
#' @export
make_sun22_hs_params <- function(par) {
  if (length(par) != 10L) {
    stop("par must have length 10: xi(2), nu(2), z(6)", call. = FALSE)
  }
  z_names <- .sun_hs_z_names(2L, 2L)
  names(par) <- c("xi1", "xi2", "nu1", "nu2", z_names)
  xi <- unname(par[c("xi1", "xi2")])
  nu <- unname(par[c("nu1", "nu2")])
  J <- .corr_from_free(unname(par[z_names]), 4L, eps = 1e-6)
  Cu <- J[1:2, 1:2, drop = FALSE]
  list(
    xi = xi,
    Omega = unname(outer(pmax(nu, 1e-12), pmax(nu, 1e-12)) * Cu),
    Delta = unname(J[1:2, 3:4, drop = FALSE]),
    tau = rep(0, 2L),
    Gamma = unname(J[3:4, 3:4, drop = FALSE])
  )
}
