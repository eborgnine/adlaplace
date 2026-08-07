#' Build SUN(4,4) parameters from a joint hyperspherical vector
#'
#' The joint correlation of \eqn{(U,V)} is formed from 28 unconstrained
#' unit-row Cholesky coordinates in row-major order.
#'
#' @param par Numeric vector of length 36: \code{xi(4)}, \code{nu(4)}, then
#'   28 \code{z} coordinates named by \code{.sun_hs_z_names(4, 4)}.
#' @return A SUN \code{dp} list with \code{xi}, \code{Omega}, \code{Delta},
#'   zero \code{tau}, and \code{Gamma}.
#' @seealso [dsun44_hs()], [make_sun44_params()], [sun44_hs_bounds()]
#' @export
make_sun44_hs_params <- function(par) {
  if (length(par) != 36L) {
    stop("par must have length 36: xi(4), nu(4), z(28)", call. = FALSE)
  }
  z_names <- .sun_hs_z_names(4L, 4L)
  names(par) <- c(paste0("xi", 1:4), paste0("nu", 1:4), z_names)
  xi <- unname(par[paste0("xi", 1:4)])
  nu <- unname(par[paste0("nu", 1:4)])
  J <- .corr_from_free(unname(par[z_names]), 8L, eps = 1e-6)
  Cu <- J[1:4, 1:4, drop = FALSE]
  list(
    xi = xi,
    Omega = unname(outer(pmax(nu, 1e-12), pmax(nu, 1e-12)) * Cu),
    Delta = unname(J[1:4, 5:8, drop = FALSE]),
    tau = rep(0, 4L),
    Gamma = unname(J[5:8, 5:8, drop = FALSE])
  )
}
