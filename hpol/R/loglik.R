#' @export
loglik <- function(
  parameters, gamma_start, 
  data, config,
  control=list()) {


  Nbeta = nrow(data$XTp)
  Ngamma = length(gamma_start)
  Sgamma0 = seq(from=Nbeta, len=Ngamma)
  Sgamma1 = Sgamma0+1



  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

  fullHessian = NULL
  if('sparsity' %in% names(config)) {
    if(any(config$sparsity$i >= Ngamma | 
      config$sparsity$j >= Ngamma)) {
      fullHessian = config$sparsity
      # this is the big hessian
      # subset to gamma-only hessian
      inSgamma = config$sparsity$i %in% Sgamma0 &
        config$sparsity$j %in% inSgamma
      config$sparsity$i = config$sparsity$i[inSgamma+1]  
      config$sparsity$j = config$sparsity$j[inSgamma+1]  
    }
  } 

  config2 = c(
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
    config = config2
  )

  result <- trustOptim::trust.optim(
    x = gamma_start,
    fn = wrappers_gamma$fn,
    gr = wrappers_gamma$gr,
    hs = wrappers_gamma$hs,
    method = "Sparse",
    control = control
  )

  result$cholHessian = Matrix::chol(result$hessian)
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
  result$parameters[-(1:Nbeta)])


parametersGamma = hpolcc::sparsityForThird(
  hessian = fullHessian, data = data)

config3 = c(
  parametersGamma,
  config2[setdiff(
      names(config2), 
      names(parametersGamma)
    )]
)
#config$verbose=TRUE

resThird = hpolcc::derivForLaplace(
  fullParameters, data, config3
  ) 

  result$resThird = resThird

  result$third = parametersGamma$full
  result$third$x = resThird$result
# TO DO: subtract off T_iik, T_jjk, H_ik, H_jk

# entry i,k is T_iik
result$deriv3diag = Matrix::sparseMatrix(
  i=parametersGamma$indexForDiag$i,
  p = parametersGamma$indexForDiag$p,
  x = resThird$forDiag, symmetric=FALSE,
  index1=FALSE, dims = rep(Nfull,2

  # to do: third derivative
#  ∂(ln(det(X))) = Tr(X^{−1} ∂X)

# gammaHat1 = gammaHat - Hinv G
    # d gammaHat1 d Theta = d Hinv / dTheta G + Hinv d G/dTheta
# dHinv = -Hinv d H/dTheta Hinv

    # d L / d theta = d L / d gammaHat *= dgammaHat d Theta


  return(result)
}
