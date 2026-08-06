#' Unit-row Cholesky (hyperspherical) factor of a correlation matrix
#'
#' Maps unconstrained free parameters to a lower-triangular factor with unit
#' Euclidean rows, so \eqn{C = B_r B_r^\top} is a correlation matrix
#' (diagonal one, no \code{cov2cor} step). Free entries use
#' \eqn{z = \tanh(\cdot)} so each partial coordinate lies in \eqn{(-1,1)}.
#'
#' @param z Numeric vector of length \code{d * (d - 1) / 2}.
#' @param d Dimension.
#' @return Lower-triangular \code{d x d} matrix with unit rows.
#' @keywords internal
.unit_row_chol <- function(z, d) {
  n_free <- as.integer(d * (d - 1L) / 2L)
  if (length(z) != n_free) {
    stop("z must have length ", n_free, " for dimension ", d, call. = FALSE)
  }
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
      Br[i, j] <- tij * sqrt(max(1 - sum_sq, 0))
      sum_sq <- sum_sq + Br[i, j]^2
    }
    Br[i, i] <- sqrt(max(1 - sum_sq, 0))
  }
  Br
}

#' Correlation matrix from unconstrained hyperspherical free parameters
#' @keywords internal
.corr_from_free <- function(z, d) {
  Br <- .unit_row_chol(z, d)
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

#' Build SUN(3,3) skew-loading matrix from diagonals, pair levels, and signed gaps
#'
#' Each off-diagonal pair shares a level \code{Lij} with a signed asymmetry
#' \code{eij}: lower triangle gets \code{Lij + eij}, upper triangle
#' \code{Lij - eij}.
#'
#' @param L11,L22,L33 Diagonal loadings.
#' @param L12,L13,L23 Pair levels for off-diagonal pairs.
#' @param e12,e13,e23 Signed lower-minus-upper gaps for those pairs.
#' @return A \eqn{3\times 3} matrix.
#' @keywords internal
.L33_from_reparam <- function(L11, L22, L33, L12, L13, L23, e12, e13, e23) {
  matrix(
    c(
      L11, L12 - e12, L13 - e13,
      L12 + e12, L22, L23 - e23,
      L13 + e13, L23 + e23, L33
    ),
    nrow = 3L,
    byrow = TRUE
  )
}

#' Invert [.L33_from_reparam()] for a given loading matrix
#' @param L A \eqn{3\times 3} matrix.
#' @return Named length-9 vector
#'   \code{c(L11,L22,L33,L12,L13,L23,e12,e13,e23)}.
#' @keywords internal
.L33_to_reparam <- function(L) {
  L <- as.matrix(L)
  if (!all(dim(L) == c(3L, 3L))) {
    stop("L must be 3 x 3", call. = FALSE)
  }
  c(
    L11 = L[1, 1], L22 = L[2, 2], L33 = L[3, 3],
    L12 = (L[2, 1] + L[1, 2]) / 2,
    L13 = (L[3, 1] + L[1, 3]) / 2,
    L23 = (L[3, 2] + L[2, 3]) / 2,
    e12 = (L[2, 1] - L[1, 2]) / 2,
    e13 = (L[3, 1] - L[1, 3]) / 2,
    e23 = (L[3, 2] - L[2, 3]) / 2
  )
}

#' Build SUN distribution parameters from a free parameter vector
#'
#' Converts the SUN(3,3) block-Cholesky parameterization used in the
#' supervision notes into the \code{dp} list expected by the \pkg{sn} package
#' (\code{xi}, \code{Omega}, \code{Delta}, \code{tau}, \code{Gamma}).
#'
#' @param par Numeric vector of length 21:
#'   \code{c(xi1, xi2, xi3, nu1, nu2, nu3, ell21, ell31, ell32,
#'   L11, L22, L33, L12, L13, L23, e12, e13, e23, a, b, c)}.
#'   Here \code{nu} are marginal standard deviations of \eqn{\Omega}
#'   (\eqn{\Omega_{jj}=\nu_j^2}), and \code{ell*}
#'   are unconstrained hyperspherical coordinates for the correlation
#'   \eqn{C_u} of the standardized latent. Skew loadings \eqn{L} act on that
#'   correlation scale (so entries in \eqn{[-1,1]} are natural); they are
#'   parameterized as three diagonals, three off-diagonal pair levels
#'   \code{L12,L13,L23}, and three signed asymmetries \code{e*} (lower
#'   \eqn{=}\code{Lij+eij}, upper \eqn{=}\code{Lij-eij}; see
#'   [.L33_from_reparam()]). \code{a}, \code{b}, \code{c} parameterize the
#'   residual correlation of \eqn{V}. Valid parameters need
#'   \eqn{\operatorname{diag}(L C_u L^\top) < 1}.
#'
#' @return A list with components \code{xi}, \code{Omega}, \code{Delta},
#'   \code{tau} (fixed at zero), and \code{Gamma}.
#'
#' @details
#' **Correlation-scale SUN, then scale.** With correlation \eqn{C_u} from the
#' unit-row Cholesky of \eqn{(\ell_{21},\ell_{31},\ell_{32})} and
#' marginal standard deviations \eqn{\nu_j>0}, writing
#' \eqn{D_u=\operatorname{diag}(\nu)},
#' \deqn{\Omega = D_u C_u D_u.}
#' Skew loadings \eqn{L} multiply \eqn{C_u} (not \eqn{\Omega}):
#' \deqn{M = L C_u L^\top,\qquad
#' \Sigma_{uv} = D_u C_u L^\top.}
#' Free \eqn{(a,b,c)} build residual correlation \eqn{C}. With
#' \eqn{D = \operatorname{diag}(1 - \operatorname{diag}(M))},
#' \deqn{R = D^{1/2} C D^{1/2},}
#' so \eqn{\Sigma_{vv} = M + R} is a correlation matrix. Then
#' \eqn{\Gamma = \Sigma_{vv}} and
#' \eqn{\Delta = D_u^{-1}\Sigma_{uv} = C_u L^\top}.
#'
#' **Skew loadings.**
#' \deqn{
#' L =
#' \begin{bmatrix}
#' L_{11} & L_{12}-e_{12} & L_{13}-e_{13} \\
#' L_{12}+e_{12} & L_{22} & L_{23}-e_{23} \\
#' L_{13}+e_{13} & L_{23}+e_{23} & L_{33}
#' \end{bmatrix}.
#' }
#'
#' @examples
#' par <- c(
#'   xi1 = 0, xi2 = 0, xi3 = 0,
#'   nu1 = 1, nu2 = 1, nu3 = 1,
#'   ell21 = 0.2, ell31 = 0.1, ell32 = 0.2,
#'   L11 = 0.4, L22 = 0.4, L33 = 0.35,
#'   L12 = 0.1, L13 = 0.05, L23 = 0.1,
#'   e12 = 0, e13 = -0.05, e23 = 0,
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
      "L11, L22, L33, L12, L13, L23, e12, e13, e23, a, b, c)"
    )
  }

  names(par) <- c(
    "xi1", "xi2", "xi3",
    "nu1", "nu2", "nu3",
    "ell21", "ell31", "ell32",
    "L11", "L22", "L33",
    "L12", "L13", "L23",
    "e12", "e13", "e23",
    "a", "b", "c"
  )

  xi <- c(par["xi1"], par["xi2"], par["xi3"])
  nu <- c(par["nu1"], par["nu2"], par["nu3"])
  Cu <- .corr_from_free(c(par["ell21"], par["ell31"], par["ell32"]), 3L)
  s_u <- pmax(as.numeric(nu), 1e-12)
  Omega <- outer(s_u, s_u) * Cu

  L <- .L33_from_reparam(
    par["L11"], par["L22"], par["L33"],
    par["L12"], par["L13"], par["L23"],
    par["e12"], par["e13"], par["e23"]
  )

  M <- L %*% Cu %*% t(L)
  R <- .residual_corr_R(M, c(par["a"], par["b"], par["c"]))

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
