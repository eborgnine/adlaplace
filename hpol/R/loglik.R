#' @export
loglik <- function(
  parameters, gamma_start, 
  data, config,
  control=list()) {


  Nbeta = nrow(data$XTp)
  Ngamma = length(gamma_start)



  beta = parameters[1:Nbeta]
  theta = parameters[-(1:Nbeta)]

  if('sparsity' %in% names(config)) {
    if(any(config$sparsity$i >= Ngamma | 
      config$sparsity$j >= Ngamma)) {
      # this is the big hessian
      # subset to gamma-only hessian
      Sgamma = seq(from=Nbeta, len=Ngamma)
      inSgamma = config$sparsity$i %in% Sgamma &
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

  result <- trust.optim(
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
    result$logdet/2 + 
    0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi

  result$gamma_hat <- result$solution

  result$parameters = parameters


  return(result)
}
