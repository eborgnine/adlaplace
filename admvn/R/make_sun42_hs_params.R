#' Build SUN(4,2) parameters from a joint hyperspherical vector
#'
#' The joint correlation of \eqn{(U,V)\in\mathbb{R}^6} is formed from 15
#' unconstrained unit-row Cholesky coordinates in row-major order.
#'
#' @param par Numeric vector of length 23: \code{xi(4)}, \code{nu(4)}, then
#'   15 \code{z} coordinates named by \code{.sun_hs_z_names(4, 2)}.
#' @return A SUN \code{dp} list with \code{xi}, \code{Omega} (\eqn{4\times 4}),
#'   \code{Delta} (\eqn{4\times 2}), zero \code{tau} (length 2), and
#'   \code{Gamma} (\eqn{2\times 2}).
#' @seealso [dsun42_hs()], [sun42_hs_bounds()]
#' @export
make_sun42_hs_params <- function(par) {
  if (length(par) != 23L) {
    stop(
      "par must have length 23: xi(4), nu(4), z(15) in order ",
      paste(.sun_hs_z_names(4L, 2L), collapse = ", "),
      call. = FALSE
    )
  }
  z_names <- .sun_hs_z_names(4L, 2L)
  names(par) <- c(paste0("xi", 1:4), paste0("nu", 1:4), z_names)
  xi <- unname(par[paste0("xi", 1:4)])
  nu <- unname(par[paste0("nu", 1:4)])
  J <- .corr_from_free(unname(par[z_names]), 6L, eps = 1e-6)
  Cu <- J[1:4, 1:4, drop = FALSE]
  list(
    xi = xi,
    Omega = unname(outer(pmax(nu, 1e-12), pmax(nu, 1e-12)) * Cu),
    Delta = unname(J[1:4, 5:6, drop = FALSE]),
    tau = rep(0, 2L),
    Gamma = unname(J[5:6, 5:6, drop = FALSE])
  )
}
