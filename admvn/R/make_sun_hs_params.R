#' Names of hyperspherical free coordinates for a joint SUN correlation
#'
#' Row-major unit-row Cholesky fill order for dimension \code{n + m}:
#' \code{z21}, \code{z31}, \code{z32}, \code{z41}, ..., \code{z\{n+m\},\{n+m-1\}}.
#'
#' @param n Dimension of \eqn{U}.
#' @param m Dimension of \eqn{V}.
#' @return Character vector of length \code{(n+m)*(n+m-1)/2}.
#' @keywords internal
.sun_hs_z_names <- function(n = 3L, m = 3L) {
  d <- as.integer(n + m)
  nm <- character(d * (d - 1L) / 2L)
  idx <- 1L
  for (i in seq_len(d)) {
    if (i == 1L) {
      next
    }
    for (j in seq_len(i - 1L)) {
      nm[[idx]] <- paste0("z", i, j)
      idx <- idx + 1L
    }
  }
  nm
}

#' Indices (1-based in the z-vector) of the independent-SN pair slots
#'
#' For square SUN(\eqn{n},\eqn{n}), the free coordinates linking
#' \eqn{V_i} to \eqn{U_i} are \eqn{z_{n+i,i}}.
#'
#' @keywords internal
.sun_hs_pair_z_idx <- function(n = 3L) {
  n <- as.integer(n)
  ## For each V-row r = n+i (i=1..n), column j=i is the i-th free entry
  ## in that row. Free count before row r is sum_{k=2}^{r-1} (k-1) = (r-1)(r-2)/2.
  vapply(seq_len(n), function(i) {
    r <- n + i
    before <- (r - 1L) * (r - 2L) / 2L
    as.integer(before + i)
  }, integer(1L))
}

#' Build SUN(3,3) dp from a joint hyperspherical parameter vector
#'
#' Alternative to [make_sun_params()]: the latent joint correlation of
#' \eqn{(U_\star, V)\in\mathbb{R}^6} is
#' \eqn{J = B B^\top} with unit-row Cholesky factor \eqn{B} from fifteen
#' unconstrained coordinates \eqn{z_{ij}} via \eqn{t_{ij}=\tanh(z_{ij})}.
#' Then \eqn{C_u=J_{UU}}, \eqn{\Delta=J_{UV}}, \eqn{\Gamma=J_{VV}}, and
#' \eqn{\Omega=D_u C_u D_u} with \eqn{D_u=\operatorname{diag}(\nu)}.
#'
#' The leading six parameters match the block-Cholesky layout
#' (\code{xi}, \code{nu}), so the C++ \code{sun_mle()} driver works
#' unchanged on hyperspherical tapes.
#'
#' @param par Numeric vector of length 21:
#'   \code{xi(3)}, \code{nu(3)}, then free \code{z} in row-major fill order
#'   \code{z21,z31,z32,z41,...,z65}. \code{nu} are marginal SDs of
#'   \eqn{\Omega} (\eqn{\Omega_{jj}=\nu_j^2}).
#'
#' @return A list with components \code{xi}, \code{Omega}, \code{Delta},
#'   \code{tau} (fixed at zero), and \code{Gamma}.
#'
#' @details
#' **Independent skew-normals.** Only the pair slots
#' \eqn{z_{4,1},z_{5,2},z_{6,3}} nonzero (equal to \eqn{\operatorname{atanh}(d_i)})
#' yields \eqn{J=\begin{bmatrix}I&D\\D&I\end{bmatrix}} with
#' \eqn{D=\operatorname{diag}(d)}. All \eqn{z=0} recovers a Gaussian.
#'
#' **Fill order.** For dimension \eqn{d=6}, free entries are filled by rows
#' \eqn{i=2,\ldots,6}, columns \eqn{j=1,\ldots,i-1}.
#'
#' @seealso [make_sun_params()], [dsun_hs()], [sun33_hs_bounds()]
#' @examples
#' par <- c(
#'   xi1 = 0, xi2 = 0, xi3 = 0,
#'   nu1 = 1, nu2 = 1, nu3 = 1,
#'   z21 = 0, z31 = 0, z32 = 0,
#'   z41 = atanh(0.4), z42 = 0, z43 = 0,
#'   z51 = 0, z52 = atanh(0.3), z53 = 0, z54 = 0,
#'   z61 = 0, z62 = 0, z63 = atanh(0.2), z64 = 0, z65 = 0
#' )
#' dp <- make_sun_hs_params(par)
#' diag(dp$Gamma)
#' @export
make_sun_hs_params <- function(par) {
  if (length(par) != 21L) {
    stop(
      "par must have length 21: xi(3), nu(3), z(15) in order ",
      paste(.sun_hs_z_names(3L, 3L), collapse = ", "),
      call. = FALSE
    )
  }
  z_names <- .sun_hs_z_names(3L, 3L)
  names(par) <- c("xi1", "xi2", "xi3", "nu1", "nu2", "nu3", z_names)

  xi <- unname(c(par["xi1"], par["xi2"], par["xi3"]))
  nu <- unname(c(par["nu1"], par["nu2"], par["nu3"]))
  z <- unname(par[z_names])

  ## eps floor keeps residual chol diagonals away from exact zero
  J <- .corr_from_free(z, 6L, eps = 1e-6)
  Cu <- J[1:3, 1:3, drop = FALSE]
  Delta <- J[1:3, 4:6, drop = FALSE]
  Gamma <- J[4:6, 4:6, drop = FALSE]

  s_u <- pmax(as.numeric(nu), 1e-12)
  Omega <- outer(s_u, s_u) * Cu

  list(
    xi = xi,
    Omega = unname(Omega),
    Delta = unname(Delta),
    tau = c(0, 0, 0),
    Gamma = unname(Gamma)
  )
}
