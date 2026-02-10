#' Reformat a Cholesky–LDLt Decomposition for Efficient Use
#'
#' Given an LDLt Cholesky factorisation returned from the C++ optimizer,
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
#'   \item{\code{D}}{A numeric vector of diagonal elements from the LDLt decomposition.}
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
#'   \item \code{halfH \%*\% H \%*\% t(halfH)}  \eqn{\approx I}.
#' }
#'
#' @param x A list with components \code{L}, \code{D}, and \code{P},
#'          typically from \code{inner_opt()} or a trust-region optimization step.
#'
#' @return A \code{dgCMatrix} giving \eqn{H^{-1/2}}.
#'
#' @examples
#' \dontrun{
#'   # Suppose inner_res$cholHessian contains the LDLt factors:
#'   halfH <- reformatChol(inner_res$cholHessian)
#'
#'   # Recover the inverse Hessian:
#'   Hinv <- Matrix::crossprod(halfH)
#'
#'   # Sanity check: Hinv %*% H approx I
#'   check <- Hinv %*% inner_res$hessian
#'   summary(as.numeric(check))
#' }
#'
#' @export
reformatChol <- function(x) {

  Linv <- Matrix::solve(x$L)
  halfDinv <- Matrix::Diagonal(ncol(x$D), (x$D@x)^(-0.5))

  # H^{-1/2} = P^T (L^{-1} D^{-1/2})
  # H^{-1/2} =  (L^{-1 T} D^{-1/2} ) P

#  halfH <- (Matrix::t(Linv) %*% halfDinv)[1 + x$P, ]

  halfH = Matrix::crossprod(Linv, halfDinv)[1 + x$P, ]
  Hinv = Matrix::tcrossprod(halfH) 

  return(list(halfH = halfH, Hinv = Hinv))
}

logLikDeriv = function(
  fullParameters,
  hessianPack,
  config, 
  adPack
) {

  derivFull = adlaplace::all_derivs(fullParameters, adPack, config)
  derivFull$hessian = do.call(Matrix::sparseMatrix, derivFull$hessian)

  Hstuff = reformatChol(hessianPack)

  Sgamma1 = seq.int(length(config$beta)+1, len=length(config$gamma))
  Sgamma0 = Sgamma1 - 1L

  dU = - Hstuff$Hinv %*% derivFull$hessian[Sgamma1, -Sgamma1]

  whichColumnsByGroup1 = lapply(
    adPack$sparsity, function(xx, refmat) {
      grad_inner_gamma = match(xx$grad_inner, Sgamma0)
      linvHere = refmat[grad_inner_gamma, ,drop=FALSE]
      which(diff(linvHere@p)>0)-1L
    }, 
    refmat = Hstuff$halfH
  )

  whichColumnsByGroup = Matrix::sparseMatrix(
    i = unlist(whichColumnsByGroup1),
    j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
    index1=FALSE,
    dims = c(length(config$gamma), length(whichColumnsByGroup1))
  )

  theTrace = adlaplace::traceHinvT(
    fullParameters, Hstuff$halfH, 
    whichColumnsByGroup,
    adPack)

  result = list(extra = list(dU = dU, trace3 = theTrace, halfH = halfH))
  result$deriv = data.frame(
    dDetUpart = as.vector(theTrace[Sgamma1] %*% dU),
    dDetTpart = theTrace[-Sgamma1])
  result$deriv$gradTheta = gradient_full[-Sgamma1]  
  result$deriv$gradU = as.vector(gradient_full[Sgamma1] %*% dU)
  result$deriv$dDet = result$deriv$dDetUpart + result$deriv$dDetTpart
  result$deriv$dL = result$deriv$dDet + result$deriv$gradU + result$deriv$gradTheta

  return(result)
}
