#' Reformat a Cholesky–LDLᵀ Decomposition for Efficient Use
#'
#' Given an LDLᵀ Cholesky factorisation returned from the C++ optimizer,
#' this function reconstructs the matrix
#' \deqn{ H^{-1/2} = P^\top L^{-1} D^{-1/2}, }
#' where the original Hessian satisfies:
#' \deqn{ H = P^\top L D L^\top P. }
#'
#' The returned \code{halfH} matrix is useful for:
#' \itemize{
#'   \item constructing \eqn{H^{-1}} via \code{crossprod(halfH)},
#'   \item whitening transformations,
#'   \item computing approximate covariance matrices in Laplace-type models.
#' }
#'
#' @details
#' The input \code{x} must be a list containing:
#' \describe{
#'   \item{\code{L}}{A lower–triangular Cholesky factor (\code{dgCMatrix}).}
#'   \item{\code{D}}{A numeric vector of diagonal elements from the LDLᵀ decomposition.}
#'   \item{\code{P}}{A permutation vector (0-based, as returned from Eigen).}
#' }
#'
#' The function computes:
#' \deqn{
#'   H^{-1}
#'   = P^\top L^{-1} D^{-1} (L^{-1})^\top P
#'   = (P^\top L^{-1} D^{-1/2}) (P^\top L^{-1} D^{-1/2})^\top.
#' }
#'
#' The object returned by \code{reformatChol()} is the matrix
#' \deqn{
#'   H^{-1/2} = P^\top L^{-1} D^{-1/2},
#' }
#' which satisfies:
#' \itemize{
#'   \item \code{crossprod(halfH)} gives \eqn{H^{-1}},
#'   \item \code{halfH \%*\% H \%*\% t(halfH)} ≈ \eqn{I}.
#' }
#'
#' @param x A list with components \code{L}, \code{D}, and \code{P},
#'          typically from \code{inner_opt()} or a trust-region optimization step.
#'
#' @return A \code{dgCMatrix} giving \(H^{-1/2}\).
#'
#' @examples
#' \dontrun{
#'   # Suppose inner_res$cholHessian contains the LDLᵀ factors:
#'   halfH <- reformatChol(inner_res$cholHessian)
#'
#'   # Recover the inverse Hessian:
#'   Hinv <- Matrix::crossprod(halfH)
#'
#'   # Sanity check: Hinv %*% H ≈ I
#'   check <- Hinv %*% inner_res$hessian
#'   summary(as.numeric(check))
#' }
#'
#' @export
reformatChol <- function(x) {

  Linv <- Matrix::solve(x$L)
  halfDinv <- Matrix::Diagonal(length(x$D), (x$D)^(-0.5))

  # H^{-1/2} = P^T (L^{-1} D^{-1/2})
  halfH <- (halfDinv %*% Linv)[, 1 + x$P]

  halfH
}