#' @export
loglik <- function(
  parameters, 
  gamma_start, 
  data, config,
  wrappers,
  control=list(), 
    deriv = c(0,1)
  ) {

  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Nparams = length(parameters) + Ngamma

  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

  if(missing(gamma_start)) gamma_start = rep(0, Ngamma)


  if(is.null(config$sparsity$third)) {
    config$sparsity = sparsityForThird(
      x=c(configInner$beta, gamma_start, configInner$theta),
      data, config)
  }


  # Optimize gamma keeping beta and theta fixed
  if(missing(wrappers)) {
        configInner = c(
      config[setdiff(
        names(config), c('beta','theta')
      )],
      list(
        beta = beta,
        theta = theta
      ))

  wrappers_gamma <- hpolcc::make_trustoptim_wrappers(
    data = data,
    config = configInner
  )
  } else {
    wrappers_gamma=wrappers
  }

  # inner opt
  result <- trustOptim::trust.optim(
    x = gamma_start,
    fn = wrappers_gamma$fn,
    gr = wrappers_gamma$gr,
    hs = wrappers_gamma$hs,
    method = "Sparse",
    control = control
  )
  result$cholHessian = Matrix::chol(result$hessian)
  invHessian = result$invHessian = Matrix::solve(result$cholHessian)

  result$logDetHessian = drop(Matrix::determinant(
    result$cholHessian, log=TRUE, sqrt=FALSE
  )$modulus)

  result$minusLogLik = result$fval +
    as.numeric(result$logDetHessian)/2 + 
    0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi

  if(all(deriv == 0)) {
    return(result$minusLogLik)
  }

  fullParameters = result$parameters = c(
    configInner$beta,
    result$solution,
    configInner$theta)

#library('hpolcc')
  resThird = derivForLaplace(
    fullParameters, data, config
  ) 

  if(FALSE) { # check hessian
    Sgamma0 = seq(from=Nbeta, len=Ngamma)
    Sgamma1 = Sgamma0+1

    bob = Matrix::Matrix(resThird$denseHessian)
    bob2 = Matrix::sparseMatrix(i=config$sparsity$second$parGamma$i,
      j = config$sparsity$second$parGamma$j, x=resThird$second, index1=FALSE)
    bob3 = Matrix::sparseMatrix(i=config$sparsity$second$parGamma$i,
      j = config$sparsity$second$parGamma$j, x=resThird$diag, index1=FALSE)

      hessianFull = Matrix::sparseMatrix(
    i = config$sparsity$second$parGamma$i,
    j = config$sparsity$second$parGamma$j,
    x = resThird$second, index1=0, dims = c(Nparams, Nparams)
  )
    range(as.matrix(hessianFull)[Sgamma1,Sgamma1] - as.matrix(result$hessian))
    range(as.matrix(bob)[Sgamma1, Sgamma1] - as.matrix(result$hessian))
    range(as.matrix(bob2)[Sgamma1, Sgamma1] - as.matrix(result$hessian))
    range((as.matrix(hessianFull) - as.matrix(resThird$denseHessian))[Sgamma1,Sgamma1])
  }

  thirdDiag = data.frame(
    i = config$sparsity$second$parGamma$j,
    j = config$sparsity$second$parGamma$j,
    k = config$sparsity$second$parGamma$i,
    x = 2*(resThird$diag - resThird$second)
  )


  thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
  thirdNonDiag$x = (
    resThird$third -
    thirdDiag$x[config$sparsity$third$ijk[,'indexKii']] -
    thirdDiag$x[config$sparsity$third$ijk[,'indexKjj']]
  )/2



  theCols = c('i','j','k', 'x')
  third = rbind(
    thirdDiag[,theCols],
    thirdNonDiag[,theCols]
  )


if(FALSE) {
  # check third deriv
  step = 0.0001;Dgamma=4
  Dorig = Dgamma + length(beta)

  gamma2 = result$solution
  gamma2[Dgamma] = result$solution[Dgamma]+step
  hessian2 = wrappers_gamma$hs(gamma2)
  D3 = (as.matrix(hessian2) - as.matrix(result$hessian))/step

  anyD = apply(third[,c('i','j','k')] == (Dorig-1), 1, any)
  autoD = third[anyD, ]

  autoD2 = as.data.frame(
    t(apply(autoD[,c('i','j','k')], 1, function(xx) sort(c(xx[xx != (Dorig-1)], rep((Dorig-1), 2))[1:2])))
  )
  autoD2$x = autoD$x


  autoDM = Matrix::sparseMatrix(i=autoD2$V1, j=autoD2$V2, x=autoD2$x, 
    symmetric=TRUE, dims = c(Nparams, Nparams), index1=FALSE)

  autoD2[autoD2$V1 == autoD2$V2, ][1:6,]
  thirdDiag[thirdDiag$k == (Dorig-1), ][1:6,]
  Matrix::diag(autoDM)[1:8]


bobDiag = Matrix::sparseMatrix(
  i=thirdDiag$i, j=thirdDiag$k, x=thirdDiag$x,
  index1=FALSE, dims = c(Nparams, Nparams)
)

  bobDiag[Dorig, Sgamma1][1:12]  

  Matrix::diag(autoDM)[1:8]
  autoDM[Sgamma1, Sgamma1][1:5,1:5]
  D3[1:5,1:5]
  # D3[1,1] is T_133


denseTkii = 2*(resThird$denseDiag - resThird$denseHessian)


range(as.matrix(bobDiag) - denseTkii)

thirdDiag[thirdDiag$k == 5,][1:8,]

(  bob = rbind(
    numer=diag(D3),
    hes=resThird$denseHessian[Dorig,Sgamma1],
    denseDiag = resThird$denseDiag[Dorig,Sgamma1],
    tdensediag=resThird$denseDiag[Sgamma1,Dorig],
    autod1=bobDiag[ Sgamma1, Dorig],
    autod2 = Matrix::diag(autoDM)[Sgamma1]
  )[,1:9])
# should be denseDiag
(bob[1,]/2 + bob[2,])[1:5]


thirdDiag[thirdDiag$i == Dorig,][1:5,]
resThird$denseDiag[1:10,Dorig]


thirdDiag[which.min(abs(thirdDiag$x - denseTkii[2,5])), ]


#  (D3/as.matrix(autoD))[1:5,1:5]
} # end testing

#  Determinant
# to do: compute in data frame.  merge and sum
  thirdList1 = mapply(function(k, third, N, Sgamma) {
    thirdHere = third[apply(third[,c('i','j','k')] == k, 1, any), ]
    newxy = t(apply(thirdHere[,c('i','j','k')], 1, function(xx) sort(c(xx[xx!=k], rep(k,2))[1:2] ) ))
   try( Matrix::sparseMatrix(i=newxy[,1], j=newxy[,2], x=thirdHere[,'x'],
      dims = rep(N, 2), symmetric=TRUE, index1=FALSE)[Sgamma, Sgamma])
  },
  k = seq(from=0, len=Nparams), MoreArgs = list(third=third, N=Nparams, Sgamma = Sgamma1)
  )


  Strace = unlist(lapply(thirdList, function(Hinv, Tuux) {
#      sum(Matrix::diag(Hinv %*% dH))
    sum((Hinv * Tuux)@x)  
  },  
  Hinv=invHessian))


  secondParGamma1 = Matrix::sparseMatrix(
    i = pmin(config$sparsity$second$parGamma$j,config$sparsity$second$parGamma$i),
    j = pmax(config$sparsity$second$parGamma$j,config$sparsity$second$parGamma$i),
    x = resThird$second,
    dims = c(Nparams, Nparams),
    index1=FALSE, symmetric=TRUE
  )
  secondParGamma = secondParGamma1[Sgamma1,-Sgamma1]

  result$extra = list(
    grad = resThird$first[-Sgamma1],
    d1 = as.vector(Strace[Sgamma1] %*% invHessian %*% secondParGamma),
    d2 = Strace[-Sgamma1]
  )

  result$gradL = result$extra$d1 + result$extra$d2 +result$extra$grad

  if(all(deriv == 1)) {
    return(result$gradL)
  } 
  result$wrappers = wrappers_gamma
  result$config = config
  result$first = resThird$first

  if(FALSE) {
    testconfig = config
    testconfig$maxDeriv=2
    hessian1 = objectiveFunctionC(
      parameters = c(beta,result$solution, theta),
      data=data, 
      config=testconfig)$hessian
    hessian1[1:6,-Sgamma1]
    secondParGamma[1:6,]
    bob = as.matrix(result$hessian) - as.matrix(hessian1[Sgamma1, Sgamma1])
  }

  return(result)
}
