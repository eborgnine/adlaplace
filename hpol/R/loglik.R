

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
    if(Nbeta == length(parameters)) {
      warning("assuming no betas, data not supplied and parameters doesn't have names")
    }
  }

  Sbeta = seq(1, len=Nbeta)
  beta = parameters[Sbeta]
  theta = parameters[setdiff(1:length(parameters), Sbeta)]
  Sgamma1 = seq(Nbeta+1, len=length(gamma_start))

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

# config$verbose=TRUE
#system.time(xx1 <- jointLogDens(gamma_start, data, config, adFun))
#system.time(xx2 <- jointLogDensDense(gamma_start, data, config, adFun))
  # inner opt
  result <- try(trustOptim::trust.optim(
    x = gamma_start,
    fn = jointLogDens,
    gr = grad,
    hs = hessian,
    method = "Sparse",
    control = control,
    data=data, config = config, adFun=adFun
  ))
  if('try-error' %in% class(result)) {
    return(list(
      minusLogLik = NA, 
      deriv=data.frame(dL = rep(NA, length(parameters)))
    ))
  }


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

  result$cholHessian = Matrix::Cholesky(result$hessian)

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

  if(identical(deriv, 0)) { # return log lik
  return(result)
}

result$invHessian = Matrix::solve(result$cholHessian)

if(missing(adFunFull)) {
  adFunFull = getAdFun(result$fullParameters, data=data, config=config)
}

result$fullHessian = hessian(
  result$fullParameters, data, config, adFunFull
)
result$fullGrad = grad(
  result$fullParameters, data, config, adFunFull
)

dU = -result$invHessian %*% result$fullHessian[Sgamma1, -Sgamma1]


cholExpand = Matrix::expand(result$cholHessian)
cholExpand$Linv = Matrix::solve(cholExpand$L)
cholExpand$LinvP = cholExpand$Linv %*% cholExpand$P
cholExpand$LinvPt = Matrix::t(cholExpand$LinvP)

linvL = as(cholExpand$LinvPt, 'lMatrix')

if(!is.null(config$first)) {
  theFirst = as(config$first, 'lMatrix')[Sgamma1,]

  whichColumnsByGroup1 = apply(theFirst, 2, function(xx) {
    linvHere = linvL[which(xx), ]
    which(diff(linvHere@p)>0)-1L
  })
  whichColumnsByGroup = Matrix::sparseMatrix(
    i = unlist(whichColumnsByGroup1),
    j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
    index1=FALSE,
    dims = dim(theFirst)
  )

  config$LinvPtColumns = whichColumnsByGroup
}

theTrace = hpolcc:::traceHinvT(
  result$fullParameters, 
  cholExpand$LinvPt, data, config, adFunFull
)

result$deriv = data.frame(
  dDetUpart = as.vector(theTrace[Sgamma1] %*% dU),
  dDetTpart = theTrace[-Sgamma1])
result$deriv$gradTheta = result$fullGrad[-Sgamma1]  
result$deriv$gradU = as.vector(result$fullGrad[Sgamma1] %*% dU)
result$deriv$dDet = result$deriv$dDetUpart + result$deriv$dDetTpart
result$deriv$dL = result$deriv$dDet + result$deriv$gradU + result$deriv$gradTheta

result$dLogLik = result$deriv$dL


return(result)
}
