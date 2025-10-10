

#' @export
loglik <- function(
  parameters, 
  gamma_start, 
  data, config,
  control=list(), 
  deriv = c(0,1)
) {

  Nbeta = nrow(data$XTp)

  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

  if(missing(gamma_start)) {
    gamma_start = rep(0, nrow(data$ATp))
  }

  if(is.null(config$sparsity$third)) {
    config$sparsity = sparsity_pattern(
      x=c(beta, gamma_start, theta),
      data, config)
  }


    configInner = c(
      config[setdiff(
        names(config), c('beta','theta')
      )],
      list(
        beta = beta,
        theta = theta
      ))


  # inner opt
  result <- trustOptim::trust.optim(
    x = gamma_start,
    fn = wrappers_gamma$fn,
    gr = wrappers_gamma$gr,
    hs = wrappers_gamma$hs,
    method = "Sparse",
    control = control,
    data=data, config = configInner
  )

  result$parameters = c(
    beta,
    theta)

  result$fullParameters = c(
    beta,
    result$solution,
    theta)

  if(identical(deriv, 0)) {

    result$cholHessian = Matrix::chol(result$hessian)
    result$invHessian = Matrix::solve(result$cholHessian)

    result$logDetHessian = drop(Matrix::determinant(
      result$cholHessian, log=TRUE, sqrt=FALSE
    )$modulus)

    result$minusLogLik = result$fval +
      as.numeric(result$logDetHessian)/2 + 
      0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi

    return(result[c('minusLogLik','solution')])
  }

  thirdRes = thirdDeriv(x=result$fullParameters, data, config)

  Sgamma1 = seq(Nbeta+1, len=length(result$solution))

  result$deriv = data.frame(
    theta = thirdRes$first[-Sgamma1],
    det = thirdRes$dDet,
    U = - as.vector(thirdRes$first[Sgamma1] %*% (thirdRes$invHessianRandom %*% thirdRes$fullHessian[Sgamma1, -Sgamma1]))
  )

  result$dLogLik =result$deriv$dL = result$deriv$theta + result$deriv$det + result$deriv$U

   if(identical(deriv, 1)) {
    return(result['deriv'])
  }


  result$extra = thirdRes
  result$config = config

  result$minusLogLik = result$fval +
    as.numeric(result$extra$logDetHessian)/2 + 
      0.5 * length(result$solution) * 1.8378770664093454835606594728  # log 2 pi



  return(result)
}
