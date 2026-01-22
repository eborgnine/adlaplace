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
  halfDinv <- Matrix::Diagonal(length(x$D), (x$D)^(-0.5))

  # H^{-1/2} = P^T (L^{-1} D^{-1/2})
  # H^{-1/2} =  (L^{-1 T} D^{-1/2} ) P
#    halfH2 <- (Matrix::t(Linv2) %*% halfDinv2) %*% cholHessian2$P1

  halfH <- (Matrix::t(Linv) %*% halfDinv)[1 + x$P, ]

  halfH
}

#' @export
logLikDeriv = function(
  x,
  inner_res,
  cholHessian=inner_res$cholHessian,
  hessian_full = inner_res$hessian_full,
  gradient_full = inner_res$gradient_full,
  config, 
  adPack,
  package = c(config$package, 'adlaplace')[1]
) {

  halfH = reformatChol(cholHessian)
  Hinv = Matrix::tcrossprod(halfH) 
  Sgamma1 = seq(length(config$beta)+1, len=nrow(Hinv))

  dU = - Hinv %*% hessian_full[Sgamma1, -Sgamma1]


  whichColumnsByGroup1 = lapply(
    config$group_inner, function(xx, refmat) {
      linvHere = refmat[1+xx$grad, ,drop=FALSE]
      which(diff(linvHere@p)>0)-1L
    }, 
    refmat = halfH
  )

  whichColumnsByGroup = Matrix::sparseMatrix(
    i = unlist(whichColumnsByGroup1),
    j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
    index1=FALSE,
    dims = c(ncol(halfH), length(whichColumnsByGroup1))
  )

  theTrace = getExportedValue(package, "traceHinvT")(
    x, halfH, 
    whichColumnsByGroup,
    config,
    adPack)



  if(FALSE){

    cholHessian2 = Matrix::expand2(Matrix::Cholesky(hessian_full[Sgamma1, Sgamma1], ldl=TRUE))
    Linv2 <- Matrix::solve(cholHessian2$L1)
    halfDinv2 <- Matrix::Diagonal(ncol(cholHessian2$D), (cholHessian2$D@x)^(-0.5))
#  halfH2 <- cholHessian2$P1. %*% (halfDinv2 %*% Linv2)
    halfH2 <- (Matrix::t(Linv2) %*% halfDinv2) %*% cholHessian2$P1

    Hinv2 = Matrix::crossprod(halfH2)

    HinvReal = Matrix::solve(hessian_full[Sgamma1, Sgamma1])

    halfH3 = Matrix::solve(Matrix::chol(hessian_full[Sgamma1, Sgamma1]))
    
    dU2 = - HinvReal %*% hessian_full[Sgamma1, -Sgamma1]


    quantile(as.matrix(Matrix::tcrossprod(halfH2) %*% hessian_full[Sgamma1, Sgamma1]) - diag(ncol(halfH)))
    quantile(as.matrix(Matrix::tcrossprod(halfH3) %*% hessian_full[Sgamma1, Sgamma1]) - diag(ncol(halfH)))
    quantile(as.matrix(Hinv %*% hessian_full[Sgamma1, Sgamma1]) - diag(ncol(halfH)))


    colFull = Matrix::sparseMatrix(
      i = rep(seq(0L, len=nrow(whichColumnsByGroup)), ncol(whichColumnsByGroup)),
      j = rep(seq(0L, len=ncol(whichColumnsByGroup)), each = nrow(whichColumnsByGroup)),
      x=1, index1=FALSE
    )
    theTrace2 = getExportedValue(package, "traceHinvT")(
      x, halfH2, 
      colFull,
      config,
      adPack)  

    theTrace3 = getExportedValue(package, "traceHinvT")(
      x, halfH3, 
      colFull,
      config,
      adPack)
    xx=cbind(as.vector(theTrace2[Sgamma1] %*% dU2), theTrace2[-Sgamma1]);cbind(xx, apply(xx, 1, sum))
    xx=cbind(as.vector(theTrace3[Sgamma1] %*% dU2), theTrace3[-Sgamma1]);cbind(xx, apply(xx, 1, sum))
    xx=cbind(as.vector(theTrace[Sgamma1] %*% dU), theTrace[-Sgamma1]);cbind(xx, apply(xx, 1, sum))




    forThird = function(x, ...) {
      adlaplace::hessian(parameters = x, ...)@x
    }
    hes1 =  numDeriv::jacobian(
      forThird,
      x = x,
      adPack = adPack, config=modifyList(config, list(verbose=FALSE))
    )


    theTraceNumeric = apply(hes1, 2, function(xx) {
      thirdHere = config$hessian_outer
      thirdHere@x = xx
      sum(Matrix::diag(Hinv2 %*% thirdHere[Sgamma1, Sgamma1]))
    })

    plot(theTraceFull, theTraceNumeric)


  }

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