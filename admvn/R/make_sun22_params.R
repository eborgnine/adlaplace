#' Build SUN(2,2) distribution parameters from a free parameter vector
#'
#' Converts the SUN(2,2) block-Cholesky parameterization into the \code{dp}
#' list expected by the \pkg{sn} package
#' (\code{xi}, \code{Omega}, \code{Delta}, \code{tau}, \code{Gamma}).
#'
#' @param par Numeric vector of length 10:
#'   \code{xi(2)}, \code{nu(2)}, \code{omega12(1)},
#'   \code{L11}, \code{L12}, \code{L21}, \code{L22}, \code{gamma12}.
#'   \code{nu} are marginal standard deviations of \eqn{\Omega}
#'   (\eqn{\Omega_{jj}=\nu_j^2}); \code{omega12} free
#'   hyperspherical coordinate for correlation \eqn{C_u}; all four entries
#'   of \eqn{L} are free on the correlation scale; \code{gamma12}
#'   parameterizes the residual correlation of \eqn{V}. Valid parameters need
#'   \eqn{\operatorname{diag}(L C_u L^\top) < 1}.
#'
#' @return A list with components \code{xi}, \code{Omega}, \code{Delta},
#'   \code{tau} (fixed at zero), and \code{Gamma}.
#'
#' @details
#' As in [make_sun_params()]: \eqn{\Omega = D_u C_u D_u} with
#' \eqn{D_u=\operatorname{diag}(\nu)},
#' \eqn{M = L C_u L^\top}, \eqn{\Sigma_{uv} = D_u C_u L^\top}, and
#' \eqn{\Sigma_{vv} = M + D^{1/2} C D^{1/2}} with unit diagonal.
#'
#' @seealso [make_sun_params()], [dsun22()], [dsun22_fun()]
#' @export
make_sun22_params <- function(par) {
  if (length(par) != 10L) {
    stop(
      "par must have length 10: xi(2), nu(2), omega12(1), ",
      "L11, L12, L21, L22, gamma12"
    )
  }
  names(par) <- c(
    "xi1", "xi2",
    "nu1", "nu2",
    "omega12",
    "L11", "L12", "L21", "L22",
    "gamma12"
  )

  xi <- unname(par[c("xi1", "xi2")])
  nu <- unname(par[c("nu1", "nu2")])
  Cu <- .corr_from_free(unname(par["omega12"]), 2L)
  s_u <- pmax(nu, 1e-12)
  Omega <- outer(s_u, s_u) * Cu

  L <- matrix(
    unname(par[c("L11", "L12", "L21", "L22")]),
    nrow = 2L,
    byrow = TRUE
  )

  M <- L %*% Cu %*% t(L)
  R <- .residual_corr_R(M, unname(par["gamma12"]))

  Delta <- Cu %*% t(L)
  Gamma <- M + R

  list(
    xi = unname(xi),
    Omega = unname(Omega),
    Delta = unname(Delta),
    tau = rep(0, 2L),
    Gamma = unname(Gamma)
  )
}
