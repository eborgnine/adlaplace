#' Build SUN(4,4) distribution parameters from a free parameter vector
#'
#' Converts the SUN(4,4) block-Cholesky parameterization into the \code{dp}
#' list expected by the \pkg{sn} package
#' (\code{xi}, \code{Omega}, \code{Delta}, \code{tau}, \code{Gamma}).
#'
#' @param par Numeric vector of length 36:
#'   \code{xi(4)}, \code{nu(4)}, \code{omega(6)}, \code{L(16)}, \code{gamma(6)}.
#'   \code{nu} are marginal standard deviations of \eqn{\Omega}
#'   (\eqn{\Omega_{jj}=\nu_j^2}); \code{omega*} free
#'   hyperspherical coordinates for correlation \eqn{C_u}; \eqn{L} acts on
#'   that correlation scale; \code{gamma*} parameterize the residual
#'   correlation of \eqn{V}. Valid parameters need
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
#' @seealso [make_sun_params()], [dsun44()], [dsun44_fun()]
#' @export
make_sun44_params <- function(par) {
  if (length(par) != 36L) {
    stop(
      "par must have length 36: xi(4), nu(4), omega(6), L(16), gamma(6)"
    )
  }
  names(par) <- c(
    "xi1", "xi2", "xi3", "xi4",
    "nu1", "nu2", "nu3", "nu4",
    "omega12", "omega13", "omega23", "omega14", "omega24", "omega34",
    "L11", "L12", "L13", "L14",
    "L21", "L22", "L23", "L24",
    "L31", "L32", "L33", "L34",
    "L41", "L42", "L43", "L44",
    "gamma12", "gamma13", "gamma23", "gamma14", "gamma24", "gamma34"
  )

  xi <- unname(par[c("xi1", "xi2", "xi3", "xi4")])
  nu <- unname(par[c("nu1", "nu2", "nu3", "nu4")])
  Cu <- .corr_from_free(
    unname(par[c(
      "omega12", "omega13", "omega23", "omega14", "omega24", "omega34"
    )]),
    4L
  )
  s_u <- pmax(nu, 1e-12)
  Omega <- outer(s_u, s_u) * Cu

  L <- matrix(
    unname(par[c(
      "L11", "L12", "L13", "L14",
      "L21", "L22", "L23", "L24",
      "L31", "L32", "L33", "L34",
      "L41", "L42", "L43", "L44"
    )]),
    nrow = 4L,
    byrow = TRUE
  )

  M <- L %*% Cu %*% t(L)
  R <- .residual_corr_R(
    M,
    unname(par[c(
      "gamma12", "gamma13", "gamma23", "gamma14", "gamma24", "gamma34"
    )])
  )

  Delta <- Cu %*% t(L)
  Gamma <- M + R

  list(
    xi = unname(xi),
    Omega = unname(Omega),
    Delta = unname(Delta),
    tau = rep(0, 4L),
    Gamma = unname(Gamma)
  )
}
