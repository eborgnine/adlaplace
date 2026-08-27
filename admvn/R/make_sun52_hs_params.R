#' Build SUN(5,2) parameters from a joint hyperspherical vector
#'
#' The joint correlation of \eqn{(U,V)\in\mathbb{R}^7} is formed from 21
#' unconstrained unit-row Cholesky coordinates in row-major order.
#'
#' @param par Numeric vector of length 31: \code{xi(5)}, \code{nu(5)}, then
#'   21 \code{z} coordinates named by \code{.sun_hs_z_names(5, 2)}.
#' @return A SUN \code{dp} list with \code{xi}, \code{Omega} (\eqn{5\times 5}),
#'   \code{Delta} (\eqn{5\times 2}), zero \code{tau} (length 2), and
#'   \code{Gamma} (\eqn{2\times 2}).
#' @seealso [dsun52_hs()], [sun52_hs_bounds()]
#' @export
make_sun52_hs_params <- function(par) {
  if (length(par) != 31L) {
    stop(
      "par must have length 31: xi(5), nu(5), z(21) in order ",
      paste(.sun_hs_z_names(5L, 2L), collapse = ", "),
      call. = FALSE
    )
  }
  z_names <- .sun_hs_z_names(5L, 2L)
  names(par) <- c(paste0("xi", 1:5), paste0("nu", 1:5), z_names)
  xi <- unname(par[paste0("xi", 1:5)])
  nu <- unname(par[paste0("nu", 1:5)])
  J <- .corr_from_free(unname(par[z_names]), 7L, eps = 1e-6)
  Cu <- J[1:5, 1:5, drop = FALSE]
  list(
    xi = xi,
    Omega = unname(outer(pmax(nu, 1e-12), pmax(nu, 1e-12)) * Cu),
    Delta = unname(J[1:5, 6:7, drop = FALSE]),
    tau = rep(0, 2L),
    Gamma = unname(J[6:7, 6:7, drop = FALSE])
  )
}
