#' @export
loglik <- function(
  parameters, 
  gamma_start, 
  data, config,
  control=list()) {

  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Sgamma0 = seq(from=Nbeta, len=Ngamma)
  Sgamma1 = Sgamma0+1

  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

  if(missing(gamma_start)) gamma_start = rep(0, Ngamma)

 if(is.null(config$sparsity$third)) {
    config$sparsity = sparsityForThird(
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

  # Optimize gamma keeping beta and theta fixed
  wrappers_gamma <- hpolcc::make_trustoptim_wrappers(
    data = data,
    config = configInner
  )

  result <- trustOptim::trust.optim(
    x = gamma_start,
    fn = wrappers_gamma$fn,
    gr = wrappers_gamma$gr,
    hs = wrappers_gamma$hs,
    method = "Sparse",
    control = control
  )

  fullParameters = c(
  beta,
  result$solution,
  theta)




# debugging
config$verbose=TRUE


resThird = derivForLaplace(
  fullParameters, data, config
  ) 




  result$cholHessian = Matrix::chol(result$hessian)
  result$invHessian = Matrix::solve(result$cholHessian)
  result$logdet = drop(Matrix::determinant(
      result$cholHessian, log=TRUE, sqrt=FALSE
    )$modulus)

  result$minusLoglik = result$fval +
    as.numeric(result$logdet)/2 + 
    0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi


  result$gamma_hat <- result$solution
  result$parameters = parameters
  Nfull = length(result$gamma_hat) + length(result$parameters)

# derivatives
if(is.null(fullHessian)) {

  fullHessian <- Matrix::Matrix(1, 
    nrow = Nfull, ncol = Nfull, 
    sparse = TRUE)
  fullHessian[
    Sgamma1, Sgamma1
  ] = result$hessian

}

fullParameters = c(
  result$parameters[1:Nbeta], 
  result$gamma_hat,
  result$parameters[-seq(1,Nbeta)])



  result$extra = resThird


  result$third = parametersGamma$full
  result$third$taylor3 = as.vector(resThird$third)




#  T i, j=param, k=gamma
#  taylor3 is T_iik + T_jjk + 2 T_ijk 

#  diag is T_iik,



#  subtract off T_iik, T_jjk,



# entry i,k of resThird$thirdis T_iik


iijIndex = 1+result$third$iInParams + nrow(resThird$diag)*result$third$paramInParams
result$third$iij = resThird$diag@x[iijIndex]
iikIndex = 1+result$third$iInParams + nrow(resThird$diag)*result$third$gammaInParams
result$third$iik = resThird$diag@x[iikIndex]
result$third$Tijk = (result$third$taylor3 - result$third$iij - result$third$iik)/2



# to do: third derivative
#  ∂(ln(det(X))) = Tr(X^{−1} ∂X)

# gammaHat1 = gammaHat - Hinv G
    # d gammaHat1 d Theta = d Hinv / dTheta G + Hinv d G/dTheta
# dHinv = -Hinv d H/dTheta Hinv
    # d gammaHat1 d Theta =  -Hinv d H/dTheta Hinv G + Hinv d G/dTheta


    # d L / d theta = d L / d gammaHat *= dgammaHat d Theta

  result$wrappers = wrappers_gamma

  return(result)
}
