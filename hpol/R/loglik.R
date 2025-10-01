

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

  Sgamma1 = seq(Nbeta+1, len=Ngamma)

  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]
  if(missing(gamma_start)) {
    gamma_start = rep(0, Ngamma)
  }

#library('hpolcc')
  if(is.null(config$sparsity$third)) {
    config$sparsity = sparsityForThird(
      x=c(beta, gamma_start, theta),
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
      # wrappers_gamma$fn(gamma_start)
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

  result$parameters = c(
    configInner$beta,
    configInner$theta)

  result$fullParameters = c(
    configInner$beta,
    result$solution,
    configInner$theta)

  if(all(deriv == 0)) {

    result$cholHessian = Matrix::chol(result$hessian)
    invHessian = result$invHessian = Matrix::solve(result$cholHessian)

    result$logDetHessian = drop(Matrix::determinant(
      result$cholHessian, log=TRUE, sqrt=FALSE
    )$modulus)

    result$minusLogLik = 
    result$fval +
    as.numeric(result$logDetHessian)/2 + 
    0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi

    return(result$minusLogLik)
  }

  thirdRes = thirdDeriv(result$fullParameters, data, config)
  thirdList = thirdRes$thirdList
  fullHessian = thirdRes$fullHessian


  cholHessianRandom = Matrix::Cholesky(fullHessian[Sgamma1, Sgamma1])
  invHessianRandom = Matrix::solve(cholHessianRandom)
  result$logDetHessian = drop(Matrix::determinant(
    cholHessianRandom, log=TRUE, sqrt=FALSE
  )$modulus)

  result$minusLogLik = result$fval +
    as.numeric(result$logDetHessian)/2 + 
      0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi


  Strace = unlist(
    lapply(thirdList, sumTrace,  
      Hinv=invHessianRandom, Sgamma1 = Sgamma1
    )
  )

  Sparams = setdiff(1:nrow(fullHessian), Sgamma1)

  StraceW = Strace[Sparams]
  StraceV = Strace[Sgamma1]
  HuuHutheta = invHessianRandom %*% fullHessian[Sgamma1, Sparams]
  VHH = Matrix::drop(StraceV %*% HuuHutheta)

  result$extra = data.frame(
    param = Sparams,
    theta = thirdRes$first[-Sgamma1],
    W = StraceW,
    VHH = VHH,
    GHH =  Matrix::drop(thirdRes$first[Sgamma1] %*% HuuHutheta) 
  )
  result$extra$grad = apply(result$extra[,c('theta','W','VHH','GHH')], 1, sum)

  if(all(deriv == 1)) {
    return(result$extra$grad)
  } 

  result$gradL = result$extra$grad
  result$wrappers = wrappers_gamma
  result$config = config
  result$first = thirdRes$first
  result$fullHessian = fullHessian
  result$dUhat = thirdRes$dUhat
  result$dH = thirdRes$dH
  result$dDet = thirdRes$dDet

  result$third = thirdList

  return(result)
}
