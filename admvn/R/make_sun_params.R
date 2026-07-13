#' Build SUN distribution parameters from a free parameter vector
#'
#' Converts the SUN(3,3) block-Cholesky parameterization used in the
#' supervision notes into the \code{dp} list expected by the \pkg{sn} package
#' (\code{xi}, \code{Omega}, \code{Delta}, \code{tau}, \code{Gamma}).
#'
#' @param par Numeric vector of length 21:
#'   \code{c(xi1, xi2, xi3, nu1, nu2, nu3, ell21, ell31, ell32,
#'   L11, ..., L33, a, b, c)}.
#'
#' @return A list with components \code{xi}, \code{Omega}, \code{Delta},
#'   \code{tau} (fixed at zero), and \code{Gamma}.
#'
#' @examples
#' par <- c(
#'   xi1 = 0, xi2 = 0, xi3 = 0,
#'   nu1 = 1, nu2 = 1, nu3 = 1,
#'   ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
#'   L11 = 1, L12 = 0.5, L13 = 1,
#'   L21 = 0.5, L22 = 1, L23 = 1.5,
#'   L31 = 0, L32 = 0.5, L33 = 1,
#'   a = 0.3, b = 0.2, c = 0.3
#' )
#' make_sun_params(par)
#'
#' @export
make_sun_params <- function(par) {
  if (length(par) != 21L) {
    stop(
      "par must have length 21: ",
      "c(xi1, xi2, xi3, nu1, nu2, nu3, ell21, ell31, ell32, ",
      "L11, L12, L13, L21, L22, L23, L31, L32, L33, a, b, c)"
    )
  }

  names(par) <- c(
    "xi1", "xi2", "xi3",
    "nu1", "nu2", "nu3",
    "ell21", "ell31", "ell32",
    "L11", "L12", "L13", "L21", "L22", "L23", "L31", "L32", "L33",
    "a", "b", "c"
  )

  xi <- c(par["xi1"], par["xi2"], par["xi3"])

  B1 <- matrix(
    c(
      1, 0, 0,
      par["ell21"], 1, 0,
      par["ell31"], par["ell32"], 1
    ),
    nrow = 3L,
    byrow = TRUE
  )
  A1 <- B1 %*% diag(c(par["nu1"], par["nu2"], par["nu3"])) %*% t(B1)

  L <- matrix(
    c(
      par["L11"], par["L12"], par["L13"],
      par["L21"], par["L22"], par["L23"],
      par["L31"], par["L32"], par["L33"]
    ),
    nrow = 3L,
    byrow = TRUE
  )

  Br <- matrix(
    c(
      1, 0, 0,
      par["a"], 1, 0,
      par["b"], par["c"], 1
    ),
    nrow = 3L,
    byrow = TRUE
  )
  R <- stats::cov2cor(Br %*% t(Br))

  Omega <- A1
  Sigma_uv <- A1 %*% t(L)
  Sigma_vv <- L %*% A1 %*% t(L) + R

  omega_u <- sqrt(diag(Omega))
  omega_v <- sqrt(diag(Sigma_vv))

  Delta <- diag(1 / omega_u) %*% Sigma_uv %*% diag(1 / omega_v)
  Gamma <- diag(1 / omega_v) %*% Sigma_vv %*% diag(1 / omega_v)

  list(
    xi = unname(xi),
    Omega = unname(Omega),
    Delta = unname(Delta),
    tau = c(0, 0, 0),
    Gamma = unname(Gamma)
  )
}
