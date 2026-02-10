

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

  result = list(extra = list(dU = dU, trace3 = theTrace, halfH = Hstuff$halfH))
  result$deriv = data.frame(
    dDetUpart = as.vector(theTrace[Sgamma1] %*% dU),
    dDetTpart = theTrace[-Sgamma1])
  result$deriv$gradTheta = derivFull$gradient[-Sgamma1]  
  result$deriv$gradU = as.vector(derivFull$gradient[Sgamma1] %*% dU)
  result$deriv$dDet = result$deriv$dDetUpart + result$deriv$dDetTpart
  result$deriv$dL = result$deriv$dDet + result$deriv$gradU + result$deriv$gradTheta

  return(result)
}
