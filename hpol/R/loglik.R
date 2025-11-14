

#' @export
loglik <- function(
  parameters, 
  gamma_start, 
  data, config,
  adFunFull,
  control=list(), 
  deriv = c(0,1),
  check=FALSE
) {

# library(hpolcc);parameters = res$parameters;gamma_start = res$gamma_start;data = res$tmb_data; config=res$config;control=res$control_inner


  if(!missing(data)) {
    Nbeta = nrow(data$XTp)
  } else {
    Nbeta = length(parameters) - sum(grepl("theta", names(parameters)))
    if(Nbeta == length(parameters)) warning("assuming no betas, data not supplied and parameters doesn't have names")
  }

  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

  if(!missing(data) & is.null(names(parameters))) {
    names(beta) = rownames(data$XTp)
    names(theta) = paste0('theta', seq(1, length(theta)))
  }

  config$beta = beta
  config$theta = theta
  
  if(missing(gamma_start)) {
    gamma_start = rep(0, nrow(data$ATp))
  }

  adFun = getAdFun(gamma_start, data=data, config=config)

  # inner opt
  result <- trustOptim::trust.optim(
    x = gamma_start,
    fn = jointLogDens,
    gr = grad,
    hs = hessian,
    method = "Sparse",
    control = control,
    data=data, config = config, adFun=adFun
  )

  if(check) {
    mleB <- try(BB::spg(par = result$solution, 
      fn = jointLogDens,
      gr = grad,
      alertConvergence=FALSE, 
      data=data, config = config, adFun=adFun,
      control = list(maxit = 1e3, M = 10, trace=TRUE))
    )
    if(!any(class(mleB) == 'try-error')) {
      result$oldSolution = result$solution
      result$solution = mleB$par
      result$hessian = hessian(mleB$par, data=data, config=config)
    }
  }

  result$parameters = c(
    beta,
    theta)

  result$fullParameters = c(
    beta,
    result$solution,
    theta)


  if(identical(deriv, 0)) { # return log lik

    # h2 = hessianDense(gamma_start, data, config)
    result$cholHessian = Matrix::Cholesky(result$hessian)
    result$invHessian = Matrix::solve(result$cholHessian)

    result$halfLogDet = try(Matrix::determinant(
      result$cholHessian, log=TRUE, sqrt=TRUE
    )$modulus)
    if(is.na(result$halfLogDet)) {
      warning("determinant of hessian is NA")
      result$halfLogDet = 1e8
    }

    result$minusLogLik = result$fval +
    as.numeric(result$halfLogDet) + 
    0.5 * length(result$solution) * 1.8378770664093454835606594728  # log 2 pi

    return(result)
  }


  if(missing(adFunFull)) {
    adFunFull = getAdFun(result$fullParameters, data=data, config=config)
  }

  thirdRes = thirdDeriv(x=result$fullParameters, data, config, adFun = adFunFull)

  Sgamma1 = seq(Nbeta+1, len=length(result$solution))


  multGammaParam = thirdRes$invHessianRandom %*% thirdRes$fullHessian[Sgamma1, -Sgamma1]

  result$deriv = data.frame(
    theta = thirdRes$first[-Sgamma1],
    det = thirdRes$dDet,
    U = - as.vector(thirdRes$first[Sgamma1] %*% multGammaParam)
  )

  result$dLogLik =result$deriv$dL = result$deriv$theta + result$deriv$det + result$deriv$U

  if(identical(deriv, 1)) {
    return(result)
  }

  result$halfLogDet = thirdRes$halfLogDet
#  result$extra = thirdRes

  result$minusLogLik = result$fval +
  as.numeric(result$halfLogDet) + 
  0.5 * length(result$solution) * 1.8378770664093454835606594728  # log 2 pi

  return(result)
}
