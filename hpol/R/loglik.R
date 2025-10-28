

#' @export
loglik <- function(
  parameters, 
  gamma_start, 
  data, config,
  control=list(), 
  deriv = c(0,1),
  check=FALSE
) {

# parameters = res$parameters;gamma_start = res$gamma_start;data = res$tmb_data; config=res$config
  Nbeta = nrow(data$XTp)

  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

    config = config[setdiff(
        names(config), c('beta','theta')
      )]

    configInner = c(
      config,
      list(
        beta = beta,
        theta = theta
      ))


  if(missing(gamma_start)) {
    gamma_start = rep(0, nrow(data$ATp))
  }


  # inner opt
  result <- trustOptim::trust.optim(
    x = gamma_start,
    fn = jointLogDens,
    gr = grad,
    hs = hessian,
    method = "Sparse",
    control = control,
    data=data, config = configInner
  )

if(check) {
mleB <- BB::spg(par = result$solution, 
    fn = jointLogDens,
    gr = grad,
  alertConvergence=FALSE, 
  data=data, config = configInner,
   control = list(maxit = 1e3, M = 10, trace=TRUE, checkGrad=TRUE))  # M = nonmonotone history
  result$solution = mleB$par
  result$hessian = hessian(mleB$par, data=data, config=configInner)
}


  if(identical(deriv, 0)) {

    result$cholHessian = Matrix::chol(result$hessian)
    result$invHessian = Matrix::solve(result$cholHessian)

    result$halfLogDet = drop(Matrix::determinant(
      result$cholHessian, log=TRUE, sqrt=TRUE
    )$modulus)

    result$minusLogLik = result$fval +
      as.numeric(result$halfLogDet) + 
      0.5 * length(result$solution) * 1.8378770664093454835606594728  # log 2 pi

    return(c(
      result[c('minusLogLik','solution', 'halfLogDet')], 
      configInner[c('beta','theta')])
    )
  }

  result$parameters = c(
    beta,
    theta)
  names(result$parameters) = c(rownames(data$XTp), paste0('theta', seq(1, length(theta))))

  result$fullParameters = c(
    beta,
    result$solution,
    theta)

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
    as.numeric(result$extra$halfLogDet) + 
      0.5 * length(result$solution) * 1.8378770664093454835606594728  # log 2 pi



  return(result)
}
