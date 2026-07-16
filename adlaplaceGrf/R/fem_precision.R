#' Assemble FEM Matern precision Q2 or Q3 from Grams
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

#' Align Gram values onto an upper-triangle CSC pattern
#' @keywords internal
align_gram_to_pattern <- function(M, p, i, n) {
  M <- methods::as(methods::as(M, "generalMatrix"), "CsparseMatrix")
  x <- numeric(length(i))
  for (col in seq_len(n) - 1L) {
    from <- p[col + 1L]
    to <- p[col + 2L] - 1L
    # Empty CSC columns have from == p[col+2]; seq.int(from, from-1) is
    # decreasing and can yield NA / invalid subscripts -- skip them.
    if (from > to) {
      next
    }
    for (pos in from:to) {
      row <- i[pos + 1L]
      x[pos + 1L] <- M[row + 1L, col + 1L]
    }
  }
  x
}

#' Build upper-triangle CSC of a structural pattern
#' @keywords internal
upper_csc_pattern <- function(S) {
  S <- methods::as(methods::as(Matrix::forceSymmetric(S), "generalMatrix"), "CsparseMatrix")
  n <- nrow(S)
  # Keep only upper triangle (row <= col) including diagonal
  U <- Matrix::triu(S)
  U <- methods::as(methods::as(U, "generalMatrix"), "CsparseMatrix")
  list(p = as.integer(U@p), i = as.integer(U@i), n = n)
}

#' Build `random_fem_ssq_*` / `random_fem_det_*` precision payload for adlaplace
#'
#' @param fem Result of [grf_bspline()] (or list with C, G, G2, optional G3).
#' @param alpha `2` or `3`.
#' @return List for `ad_data@precision`: Grams, chol pattern, and Q CSC
#'   coefficients aligned for on-tape assembly.
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
  pat <- upper_csc_pattern(struct)
  out <- list(
    C = C,
    G = G,
    G2 = G2,
    chol = chol,
    alpha = alpha,
    Q_p = pat$p,
    Q_i = pat$i,
    C_x = align_gram_to_pattern(C, pat$p, pat$i, pat$n),
    G_x = align_gram_to_pattern(G, pat$p, pat$i, pat$n),
    G2_x = align_gram_to_pattern(G2, pat$p, pat$i, pat$n)
  )
  if (!is.null(G3)) {
    out$G3 <- G3
    out$G3_x = align_gram_to_pattern(G3, pat$p, pat$i, pat$n)
  } else {
    out$G3_x <- rep(0, length(pat$i))
  }
  out
}
