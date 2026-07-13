#' Symbolic LDL pattern for a sparse FEM precision structure
#'
#' Builds the permutation and unit-lower pattern required by
#' `adlaplace`'s `chol_update` for a precision with the same sparsity as
#' `Q_structure` (typically the union of FEM Gram patterns).
#'
#' @param Q_structure Symmetric sparse matrix with the structural nonzeros of
#'   `Q(θ)` (numeric values only need to make the matrix positive definite for
#'   symbolic analysis; diagonals are inflated if needed).
#' @return A list with `perm` (0-based), `L1` (`dtCMatrix` / unit lower), and
#'   `perm_inv` (0-based), suitable for the `chol` slot of a `random_fem_*`
#'   precision payload.
#' @export
fem_chol_pattern <- function(Q_structure) {
  Q <- methods::as(methods::as(Q_structure, "generalMatrix"), "CsparseMatrix")
  n <- nrow(Q)
  if (n != ncol(Q)) {
    stop("Q_structure must be square")
  }
  # Ensure a usable PD matrix for Cholesky analysis
  Q <- Matrix::forceSymmetric(Q)
  d <- Matrix::diag(Q)
  if (any(d <= 0) || any(!is.finite(d))) {
    Matrix::diag(Q) <- pmax(abs(d), 1) + 1
  } else {
    Matrix::diag(Q) <- d + 1
  }
  ch <- Matrix::Cholesky(Q, perm = TRUE, LDL = TRUE)
  L <- Matrix::expand2(ch)$L1
  perm0 <- as.integer(ch@perm)
  list(
    perm = perm0,
    perm_inv = as.integer(order(ch@perm) - 1L),
    L1 = L
  )
}

#' Structural nonzero pattern for Q2 or Q3 from Grams
#' @export
fem_Q_structure <- function(C, G, G2, G3 = NULL) {
  S <- Matrix::drop0(abs(C) + abs(G) + abs(G2))
  if (!is.null(G3)) {
    S <- Matrix::drop0(S + abs(G3))
  }
  methods::as(methods::as(Matrix::forceSymmetric(S), "generalMatrix"), "CsparseMatrix")
}
