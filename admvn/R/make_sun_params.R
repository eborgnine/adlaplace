#' Unit-row Cholesky (hyperspherical) factor of a correlation matrix
#'
#' Maps unconstrained free parameters to a lower-triangular factor with unit
#' Euclidean rows, so \eqn{C = B_r B_r^\top} is a correlation matrix
#' (diagonal one, no \code{cov2cor} step). Free entries use
#' \eqn{z = \tanh(\cdot)} so each partial coordinate lies in \eqn{(-1,1)}.
#'
#' @param z Numeric vector of length \code{d * (d - 1) / 2}.
#' @param d Dimension.
#' @param eps Floor for residual row-norm squares (default \code{0}).
#' @return Lower-triangular \code{d x d} matrix with unit rows.
#' @keywords internal
.unit_row_chol <- function(z, d, eps = 0) {
  n_free <- as.integer(d * (d - 1L) / 2L)
  if (length(z) != n_free) {
    stop("z must have length ", n_free, " for dimension ", d, call. = FALSE)
  }
  eps <- max(as.numeric(eps), 0)
  Br <- diag(d)
  idx <- 1L
  for (i in seq_len(d)) {
    if (i == 1L) {
      next
    }
    sum_sq <- 0
    for (j in seq_len(i - 1L)) {
      tij <- tanh(z[[idx]])
      idx <- idx + 1L
      Br[i, j] <- tij * sqrt(max(1 - sum_sq, eps))
      sum_sq <- sum_sq + Br[i, j]^2
    }
    Br[i, i] <- sqrt(max(1 - sum_sq, eps))
  }
  Br
}

#' Correlation matrix from unconstrained hyperspherical free parameters
#' @keywords internal
.corr_from_free <- function(z, d, eps = 0) {
  Br <- .unit_row_chol(z, d, eps = eps)
  Br %*% t(Br)
}

#' Covariance from standard deviations and unconstrained correlation free params
#'
#' \eqn{A = D C D} with \eqn{D = \operatorname{diag}(\nu)} (\code{nu} are
#' standard deviations) and \eqn{C} from [.corr_from_free()].
#' @keywords internal
.cov_from_sd_corr <- function(nu, z) {
  nu <- as.numeric(nu)
  d <- length(nu)
  C <- .corr_from_free(z, d)
  s <- pmax(nu, 1e-12)
  outer(s, s) * C
}

#' @keywords internal
.cov_from_var_corr <- function(nu, z) {
  .cov_from_sd_corr(sqrt(pmax(as.numeric(nu), 1e-12)), z)
}

#' Invert unit-row Cholesky + tanh map: correlation -> free parameters
#' @keywords internal
.corr_free_from_C <- function(C, eps = 1e-8) {
  C <- as.matrix(C)
  d <- nrow(C)
  ## Lower Cholesky of a correlation has unit rows
  L <- t(chol((C + t(C)) / 2))
  z <- numeric(d * (d - 1L) / 2L)
  idx <- 1L
  for (i in seq_len(d)) {
    if (i == 1L) {
      next
    }
    sum_sq <- 0
    for (j in seq_len(i - 1L)) {
      rem <- max(1 - sum_sq, eps)
      tij <- L[i, j] / sqrt(rem)
      tij <- min(max(tij, -1 + eps), 1 - eps)
      z[idx] <- atanh(tij)
      idx <- idx + 1L
      sum_sq <- sum_sq + L[i, j]^2
    }
  }
  z
}

#' Residual covariance so that M + R has unit diagonal
#'
#' Given \eqn{M = L C_u L^\top} and unconstrained hyperspherical coordinates
#' \code{z} for a residual correlation \eqn{C}, returns
#' \eqn{R = D^{1/2} C D^{1/2}} with \eqn{D = \operatorname{diag}(1 - \operatorname{diag}(M))}.
#' Requires \eqn{\operatorname{diag}(M) < 1} for a valid residual variance
#' (entries are clamped at a small positive floor for numerical safety).
#'
#' @param M Symmetric \eqn{d \times d} matrix.
#' @param z Free residual-correlation parameters, length \code{d*(d-1)/2}.
#' @param eps Floor for residual variances (default \code{1e-12}).
#' @return Symmetric residual covariance \eqn{R}.
#' @keywords internal
.residual_corr_R <- function(M, z, eps = 1e-12) {
  M <- as.matrix(M)
  d <- nrow(M)
  C <- .corr_from_free(z, d)
  rem <- pmax(1 - diag(M), eps)
  s <- sqrt(rem)
  outer(s, s) * C
}

#' Build SUN distribution parameters from a free parameter vector
#'
#' Converts the SUN(3,3) block-Cholesky parameterization used in the
#' supervision notes into the \code{dp} list expected by the \pkg{sn} package
#' (\code{xi}, \code{Omega}, \code{Delta}, \code{tau}, \code{Gamma}).
#'
#' @param par Numeric vector of length 21:
#'   \code{c(xi1, xi2, xi3, nu1, nu2, nu3, omega12, omega13, omega23,
#'   L11, L12, L13, L21, L22, L23, L31, L32, L33, gamma12, gamma13, gamma23)}.
#'   Here \code{nu} are marginal standard deviations of \eqn{\Omega}
#'   (\eqn{\Omega_{jj}=\nu_j^2}), and \code{omega*}
#'   are unconstrained hyperspherical coordinates for the correlation
#'   \eqn{C_u} of the standardized latent. Skew loadings \eqn{L} act on that
#'   correlation scale (all nine entries free; diagonals typically in
#'   \eqn{[-1,1]}, off-diagonals in \eqn{[-0.5,0.5]}). \code{gamma*}
#'   parameterize the residual correlation of \eqn{V}. Valid
#'   parameters need \eqn{\operatorname{diag}(L C_u L^\top) < 1}.
#'
#' @return A list with components \code{xi}, \code{Omega}, \code{Delta},
#'   \code{tau} (fixed at zero), and \code{Gamma}.
#'
#' @details
#' **Correlation-scale SUN, then scale.** With correlation \eqn{C_u} from the
#' unit-row Cholesky of \eqn{(\omega_{12},\omega_{13},\omega_{23})} and
#' marginal standard deviations \eqn{\nu_j>0}, writing
#' \eqn{D_u=\operatorname{diag}(\nu)},
#' \deqn{\Omega = D_u C_u D_u.}
#' Skew loadings \eqn{L} multiply \eqn{C_u} (not \eqn{\Omega}):
#' \deqn{M = L C_u L^\top,\qquad
#' \Sigma_{uv} = D_u C_u L^\top.}
#' Free \eqn{(\gamma_{12},\gamma_{13},\gamma_{23})} build residual
#' correlation \eqn{C}. With
#' \eqn{D = \operatorname{diag}(1 - \operatorname{diag}(M))},
#' \deqn{R = D^{1/2} C D^{1/2},}
#' so \eqn{\Sigma_{vv} = M + R} is a correlation matrix. Then
#' \eqn{\Gamma = \Sigma_{vv}} and
#' \eqn{\Delta = D_u^{-1}\Sigma_{uv} = C_u L^\top}.
#'
#' **Skew loadings.** All nine entries of \eqn{L} are free:
#' \deqn{
#' L =
#' \begin{bmatrix}
#' L_{11} & L_{12} & L_{13} \\
#' L_{21} & L_{22} & L_{23} \\
#' L_{31} & L_{32} & L_{33}
#' \end{bmatrix}.
#' }
#'
#' @examples
#' par <- c(
#'   xi1 = 0, xi2 = 0, xi3 = 0,
#'   nu1 = 1, nu2 = 1, nu3 = 1,
#'   omega12 = 0.2, omega13 = 0.1, omega23 = 0.2,
#'   L11 = 0.4, L12 = 0.1, L13 = 0.1,
#'   L21 = 0.1, L22 = 0.4, L23 = 0.1,
#'   L31 = 0, L32 = 0.1, L33 = 0.35,
#'   gamma12 = 0.3, gamma13 = 0.2, gamma23 = 0.3
#' )
#' make_sun_params(par)
#'
#' @export
make_sun_params <- function(par) {
  if (length(par) != 21L) {
    stop(
      "par must have length 21: ",
      "c(xi1, xi2, xi3, nu1, nu2, nu3, omega12, omega13, omega23, ",
      "L11, L12, L13, L21, L22, L23, L31, L32, L33, ",
      "gamma12, gamma13, gamma23)"
    )
  }

  names(par) <- c(
    "xi1", "xi2", "xi3",
    "nu1", "nu2", "nu3",
    "omega12", "omega13", "omega23",
    "L11", "L12", "L13", "L21", "L22", "L23", "L31", "L32", "L33",
    "gamma12", "gamma13", "gamma23"
  )

  xi <- c(par["xi1"], par["xi2"], par["xi3"])
  nu <- c(par["nu1"], par["nu2"], par["nu3"])
  Cu <- .corr_from_free(c(par["omega12"], par["omega13"], par["omega23"]), 3L)
  s_u <- pmax(as.numeric(nu), 1e-12)
  Omega <- outer(s_u, s_u) * Cu

  L <- matrix(
    c(
      par["L11"], par["L12"], par["L13"],
      par["L21"], par["L22"], par["L23"],
      par["L31"], par["L32"], par["L33"]
    ),
    nrow = 3L,
    byrow = TRUE
  )

  M <- L %*% Cu %*% t(L)
  R <- .residual_corr_R(M, c(par["gamma12"], par["gamma13"], par["gamma23"]))

  Delta <- Cu %*% t(L)
  Gamma <- M + R

  list(
    xi = unname(xi),
    Omega = unname(Omega),
    Delta = unname(Delta),
    tau = c(0, 0, 0),
    Gamma = unname(Gamma)
  )
}
