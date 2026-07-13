#' Assemble FEM Matérn precision Q2 or Q3 from Grams
#'
#' @param kappa,tau Positive SPDE parameters.
#' @param C,G,G2,G3 Sparse Grams from [grf_bspline()].
#' @param alpha `2` or `3`.
#' @return Sparse `dgCMatrix` precision.
#' @export
fem_precision <- function(kappa, tau, C, G, G2, G3 = NULL, alpha = 2L) {
  kappa <- as.numeric(kappa)
  tau <- as.numeric(tau)
  if (kappa <= 0 || tau <= 0) {
    stop("kappa and tau must be positive")
  }
  alpha <- as.integer(alpha)
  if (alpha == 2L) {
    Q <- tau^2 * (kappa^4 * C + 2 * kappa^2 * G + G2)
  } else if (alpha == 3L) {
    if (is.null(G3)) {
      stop("G3 is required for alpha = 3")
    }
    Q <- tau^2 * (kappa^6 * C + 3 * kappa^4 * G + 3 * kappa^2 * G2 + G3)
  } else {
    stop("alpha must be 2 or 3")
  }
  methods::as(
    methods::as(Matrix::drop0(Matrix::forceSymmetric(Q)), "generalMatrix"),
    "CsparseMatrix"
  )
}

#' Build `random_fem_2` / `random_fem_3` precision payload for adlaplace
#'
#' @param fem Result of [grf_bspline()] (or list with C, G, G2, optional G3).
#' @param alpha `2` or `3`.
#' @return List for `ad_data@precision`: Grams plus `chol` pattern.
#' @export
fem_precision_payload <- function(fem, alpha = 2L) {
  alpha <- as.integer(alpha)
  C <- fem$C
  G <- fem$G
  G2 <- fem$G2
  G3 <- if (alpha >= 3L) fem$G3 else NULL
  if (alpha >= 3L && is.null(G3)) {
    stop("fem$G3 is required for alpha = 3 (use degree >= 3 in grf_bspline)")
  }
  struct <- fem_Q_structure(C, G, G2, G3)
  chol <- fem_chol_pattern(struct)
  out <- list(C = C, G = G, G2 = G2, chol = chol, alpha = alpha)
  if (!is.null(G3)) {
    out$G3 <- G3
  }
  out
}
