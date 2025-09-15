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
  Ntheta = length(theta)
  Nparams = Nbeta + Ngamma + Ntheta

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

resThird = derivForLaplace(
  fullParameters, data, config
  ) 

thirdDiag = data.frame(
  i = config$sparsity$second$parGamma$i,
  k = config$sparsity$second$parGamma$j,
  Tkii = resThird$diag
)
thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
thirdNonDiag$taylor = resThird$third
thirdNonDiag = merge(thirdNonDiag, thirdDiag, all.x=TRUE, all.y=FALSE)
names(thirdDiag) = gsub("i", "j", names(thirdDiag))
thirdNonDiag = merge(thirdNonDiag, thirdDiag, all.x=TRUE, all.y=FALSE)
# taylor is   //  T_iik + T_jjk + 2 T_ijk 
thirdNonDiag$x = (thirdNonDiag$taylor - 
  thirdNonDiag$Tkii - thirdNonDiag$Tkjj)/2

thirdDiag$i = thirdDiag$j
names(thirdDiag) = gsub("Tk..", "x", names(thirdDiag))

theCols = c('i','j','k', 'x')
third = rbind(
  thirdDiag[,theCols],
  thirdNonDiag[,theCols]
)
thirdList1 = split(third, third$k)
thirdList = lapply(thirdList1, function(xx, dims, Sgamma) {
  Matrix::sparseMatrix(i=xx$i, j=xx$j, x=xx$x, symmetric=TRUE, dims=dims, index1=FALSE)[Sgamma, Sgamma]
}, dims=c(Nparams,Nparams), Sgamma=Sgamma1)

secondParGamma1 = Matrix::sparseMatrix(
  i = pmin(config$sparsity$second$parGamma$j,config$sparsity$second$parGamma$i),
  j = pmax(config$sparsity$second$parGamma$j,config$sparsity$second$parGamma$i),
  x = resThird$second,
  dims = c(Nparams, Nparams),
  index1=FALSE, symmetric=TRUE
)
secondParGamma = secondParGamma1[Sgamma1,-Sgamma1]

cholHessian = Matrix::chol(result$hessian)
invHessian = Matrix::solve(cholHessian)

#  Determinant
# to do: compute in data frame.  merge and sum
dDetW = unlist(lapply(thirdList, function(dH, Hinv, HinvGamma) {
#      sum(Matrix::diag(Hinv %*% dH))
      sum((Hinv * cH)@x)  
  },  
  Hinv=invHessian, 
  HinvGamma = drop(invHessian %*% result$solution)))

# to do: V trace part

dGammapart = invHessian %*% secondParGamma

# gammaHat1 = gammaHat - Hinv G
    # d gammaHat1 d Theta = d Hinv / dTheta G + Hinv d G/dTheta
# dHinv = -Hinv d H/dTheta Hinv
    # d gammaHat1 d Theta =  -Hinv d H/dTheta Hinv G + Hinv d G/dTheta


    # d L / d theta = d L / d gammaHat *= dgammaHat d Theta


# to do: compute full hessian, check result$hessian and secondParGamma
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




  result$wrappers = wrappers_gamma

  return(result)
}
