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

  if(FALSE) { # check
    Sgamma0 = seq(from=Nbeta, len=Ngamma)
    Sgamma1 = Sgamma0+1

   # bob = Matrix::Matrix(resThird$denseHessian)
    bob2 = Matrix::sparseMatrix(i=config$sparsity$second$parGamma$i,
      j = config$sparsity$second$parGamma$j, x=resThird$second, index1=FALSE)
    bob3 = Matrix::sparseMatrix(i=config$sparsity$second$parGamma$i,
      j = config$sparsity$second$parGamma$j, x=resThird$diag, index1=FALSE)
   # range(as.matrix(bob)[Sgamma1, Sgamma1] - as.matrix(result$hessian))
    range(as.matrix(bob2)[Sgamma1, Sgamma1] - as.matrix(result$hessian))
  }


  thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
  thirdNonDiag$x = (
    resThird$third -
    resThird$diag[config$sparsity$third$ijk[,'indexKii']] -
    resThird$diag[config$sparsity$third$ijk[,'indexKjj']]
  )/2

  thirdDiag = data.frame(
    i = config$sparsity$second$parGamma$i,
    j = config$sparsity$second$parGamma$i,
    k = config$sparsity$second$parGamma$j,
    x = resThird$diag
  )

  theCols = c('i','j','k', 'x')
  third = rbind(
    thirdDiag[,theCols],
    thirdNonDiag[,theCols]
  )

#  Determinant
# to do: compute in data frame.  merge and sum
  thirdList1 = split(third, third$k)
  thirdList = lapply(thirdList1, function(xx, dims, Sgamma) {
    Matrix::sparseMatrix(i=xx$i, j=xx$j, x=xx$x, symmetric=TRUE, dims=dims, index1=FALSE)[Sgamma, Sgamma]
  }, dims=c(Nparams,Nparams), Sgamma=Sgamma1)

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
