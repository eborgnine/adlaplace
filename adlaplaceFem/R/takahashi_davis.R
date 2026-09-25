#' Sparse inverse subset of a symmetric positive-definite matrix
#'
#' Entries of \eqn{Q^{-1}} on the sparsity pattern of the Cholesky factor, in
#' the original variable order. Positions the factor does not store are absent,
#' not zeros of the inverse.
#'
#' C++ callers `LinkingTo: adlaplaceFem` and include
#' `adlaplaceFem/takahashi_davis.hpp`. That function takes the unit-lower
#' factor \eqn{L} and the LDL diagonal \eqn{D}. A caller that has only
#' \eqn{Q} factors it with `chol_update_csc` from adlaplace first.
#'
#' @param Q Symmetric sparse positive-definite matrix.
#' @param chol Optional list from [fem_chol_pattern()]. When `NULL`, the
#'   symbolic factor is computed from `Q`.
#' @return A `dgCMatrix` of the selected inverse.
#' @export
takahashi_davis <- function(Q, chol = NULL) {
  Q <- methods::as(Q, "CsparseMatrix")
  if (nrow(Q) != ncol(Q)) {
    stop("takahashi_davis: Q must be square")
  }
  Q <- Matrix::forceSymmetric(Q)
  if (is.null(chol)) {
    chol <- fem_chol_pattern(Q)
  }
  if (is.null(chol$perm) || is.null(chol$L1)) {
    stop("chol must be the list returned by fem_chol_pattern()")
  }
  U <- Matrix::triu(Q)
  U <- methods::as(methods::as(U, "generalMatrix"), "CsparseMatrix")
  takahashi_davis_matrix(U, chol$perm, chol$L1)
}
