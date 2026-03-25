

reformatChol <- function(x) {

  Linv <- Matrix::solve(x$L1)
  halfDinv <- Matrix::Diagonal(length(x$D@x), (x$D@x)^(-0.5))

  # H^{-1/2} = P^T (L^{-1} D^{-1/2})
  # H^{-1/2} =  (L^{-1 T} D^{-1/2} ) P

#  halfH <- (Matrix::t(Linv) %*% halfDinv)[1 + x$P, ]

  halfH = Matrix::crossprod(Linv, halfDinv)[x$P1@perm, ]
  Hinv = Matrix::tcrossprod(halfH) 

  return(list(halfH = halfH, Hinv = Hinv))
}

logLikDeriv = function(
  full_parameters,
  hessianPack,
  grad,
  config, 
  adFun
) {


  Hstuff = reformatChol(hessianPack$cholInner)

  Sgamma1 = seq.int(length(config$beta)+1, len=length(config$gamma))
  Sgamma0 = Sgamma1 - 1L
  # tocheck = Hstuff$Hinv %*% hessianPack$H[Sgamma1, Sgamma1]


  # to do: grad_inner_gamma computed when ADfun is created
  whichColumnsByGroup1 = lapply(
    adFun$sparsity, function(xx, refmat) {
      grad_inner_gamma = match(xx$grad_inner, Sgamma0)
      linvHere = refmat[grad_inner_gamma, ,drop=FALSE]
      which(diff(linvHere@p)>0)-1L # which columns have at least one non-zero
    }, 
    refmat = Hstuff$halfH
  )

  whichColumnsByGroup = Matrix::sparseMatrix(
    i = unlist(whichColumnsByGroup1),
    j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
    index1=FALSE,
    dims = c(length(config$gamma), length(whichColumnsByGroup1))
  )

  # need to pass num threads
  theTrace = adlaplace::traceHinvT(
    full_parameters, Hstuff$halfH, 
    whichColumnsByGroup,
    adFun, 
    c(config$num_threads, 1L)[1]
  )

  dU = - Hstuff$Hinv %*% hessianPack$H[Sgamma1, -Sgamma1]

  result = list(extra = list(dU = dU, trace3 = theTrace, halfHinv = Hstuff$halfH))


  # until now quantities are deriviative of negative log likelihood.
  # the inner_opt functions produce derivatives for neg log dens
  # below are derivatives for log likelkihood
  result$deriv = data.frame(
    dDetUpart = -as.vector(theTrace[Sgamma1] %*% dU),
    dDetTpart = -theTrace[-Sgamma1])
  result$deriv$gradTheta = -grad[-Sgamma1]  
  result$deriv$gradU = as.vector(-grad[Sgamma1] %*% dU)
  result$deriv$dDet = result$deriv$dDetUpart + result$deriv$dDetTpart
  result$deriv$dL = result$deriv$gradTheta - result$deriv$dDet # + result$deriv$gradU 

  return(result)
}
