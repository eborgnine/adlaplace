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

    result$parameters = c(
      configInner$beta,
      configInner$theta)

    fullParameters = c(
      configInner$beta,
      result$solution,
      configInner$theta)

    resThird = derivForLaplace(
      fullParameters, data, config
    ) 

    # T_iik
    thirdDiag = data.frame(
      i = config$sparsity$second$nonSymmetric$j,
      j = config$sparsity$second$nonSymmetric$j,
      k = config$sparsity$second$nonSymmetric$i,
      x = 2*(resThird$diag - resThird$second)
    )
    fullHessian = Matrix::forceSymmetric(Matrix::sparseMatrix(
      j = config$sparsity$second$nonSymmetric$j,
      i = config$sparsity$second$nonSymmetric$i,
      x = resThird$second,
      dims = rep(length(resThird$first), 2), index1=FALSE
    ))

  # taylor is   T_iik/2 + T_jjk/2 + T_ijk 
    thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
    thirdNonDiag$x = (
      resThird$third - #thirdDiag$x[Mik] - thirdDiag$x[Mjk]
      (thirdDiag$x[config$sparsity$third$ijk[,'diagIk1']] +
      thirdDiag$x[config$sparsity$third$ijk[,'diagJk1']])/2
    )



    theCols = c('i','j','k', 'x')
    third = rbind(
      thirdDiag[,theCols],
      thirdNonDiag[,theCols]
    )

    thirdList = mapply(
      function(k, third, N, Sgamma) {
        thirdHere = third[apply(third[,c('i','j','k')] == k, 1, any), ]
        if(!nrow(thirdHere)) return(NULL)
          newxy = t(apply(thirdHere[,c('i','j','k')], 1, function(xx) sort(c(xx[xx!=k], rep(k,2))[1:2] ) ))
        try(Matrix::sparseMatrix(i=newxy[,1], j=newxy[,2], x=thirdHere[,'x'],
          dims = rep(N, 2), symmetric=TRUE, index1=FALSE)) #[Sgamma, Sgamma]
      },
      k = seq(from=0, len=length(config$sparsity$full$p)), 
      MoreArgs = list(third=third, N=length(config$sparsity$full$p), Sgamma = Sgamma1)
    )


    if(FALSE) { # check hessian
    Sgamma0 = seq(from=Nbeta, len=Ngamma)
    Sgamma1 = Sgamma0+1



    bob = Matrix::Matrix(resThird$denseHessian)
    bob2 = Matrix::sparseMatrix(i=config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$j, x=resThird$second, 
      dims = dim(resThird$denseHessian), index1=FALSE)
    bob3 = Matrix::sparseMatrix(i=config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$j, x=resThird$diag, dims = dim(resThird$denseHessian), index1=FALSE)

    hessianFull = Matrix::sparseMatrix(
      i = config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$j,
      x = resThird$second, index1=0, dims = c(Nparams, Nparams)
    )
    hessianFull2=  as.matrix(Matrix::forceSymmetric(hessianFull, "L"))
    quantile(hessianFull2[Sgamma1,Sgamma1] - as.matrix(result$hessian))
    quantile(as.matrix(bob)[Sgamma1, Sgamma1] - as.matrix(result$hessian))
    quantile(as.matrix(Matrix::forceSymmetric(bob2, "L"))[Sgamma1, Sgamma1] - as.matrix(result$hessian))
    quantile((as.matrix(Matrix::forceSymmetric(hessianFull, "L")) - as.matrix(resThird$denseHessian))[Sgamma1,Sgamma1])

    quantile(as.matrix(bob2) - resThird$denseHessian)
    quantile(as.matrix(result$hessian) - resThird$denseHessian[Sgamma1, Sgamma1])

  }


  if(FALSE) {



    denseThirdIndex = Matrix::sparseMatrix(
      i = config$sparsity$third$ijk$i, 
      j=config$sparsity$third$ijk$j,
      x=(config$sparsity$third$ijk$k),
      index1=FALSE)

    options(width=180,digits=3)
    resThird$denseDiag[1:11,1:11]
    resThird$denseHessian[1:11,1:11]

    denseTkii = 2*(resThird$denseDiag - resThird$denseHessian)
    denseTkii[1:11,1:11]

    bobTH = Matrix::sparseMatrix(
      i = config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$j,
      x = resThird$diag , dims = c(Nparams, Nparams), index1=0, symmetric=FALSE
    )
    bobH = Matrix::sparseMatrix(
      i = config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$j,
      x = resThird$second , dims = c(Nparams, Nparams), index1=0, symmetric=FALSE
    )
    bobT = 2*(bobTH - bobH)

    quantile(as.matrix(bobTH) - resThird$denseDiag)

  # check third deriv
    step = 0.0001;Dgamma=1
    Dorig = Dgamma + length(beta)

    gamma1 = gamma2 = result$solution
    gamma2[Dgamma] = result$solution[Dgamma]+step
    hessianA = wrappers_gamma$hs(gamma1)
    hessianB = wrappers_gamma$hs(gamma2)

    (Matrix::determinant(hessianB)$modulus - Matrix::determinant(hessianA)$modulus)/step

    D3 = (hessianB - hessianA)/step
    trace(solve(hessianA, ))

    Sindex = 1:10

    thirdList[[Dorig]][2+Sindex,2+Sindex]
    D3[Sindex,Sindex]


    anyD = apply(third[,c('i','j','k')] == (Dorig-1), 1, any)
    autoD = third[anyD, ]

    autoD2 = as.data.frame(
      t(apply(autoD[,c('i','j','k')], 1, function(xx) sort(c(xx[xx != (Dorig-1)], rep((Dorig-1), 2))[1:2])))
    )
    autoD2$x = autoD$x

    autoDM = Matrix::sparseMatrix(i=autoD2$V1, j=autoD2$V2, x=autoD2$x, 
      symmetric=TRUE, dims = c(Nparams, Nparams), index1=FALSE)[Sgamma1,Sgamma1]

    autoD2[autoD2$V1 == autoD2$V2, ][1:6,]
    thirdDiag[thirdDiag$k == (Dorig-1), ][1:6,]
    Matrix::diag(autoDM)[1:8]


    bobDiag = Matrix::sparseMatrix(
      i=thirdDiag$i, j=thirdDiag$k, x=thirdDiag$x,
      index1=FALSE, dims = c(Nparams, Nparams)
    )

    bobDiag[Dorig, Sgamma1][1:12]  

    Matrix::diag(autoDM)[1:8]
    autoDM[Sgamma1, Sgamma1][1:5,1:5]
    D3[1:5,1:5]
  # D3[1,1] is T_133





    range(as.matrix(bobDiag) - denseTkii)

    thirdDiag[thirdDiag$k == 5,][1:8,]

# checking diagonals
    (  bob = rbind(
      numer=diag(D3),
      hes=resThird$denseHessian[Dorig,Sgamma1],
      denseDiag = resThird$denseDiag[Dorig,Sgamma1],
      tdensediag=denseTkii[Dorig,Sgamma1],
      autod1=bobDiag[ Sgamma1, Dorig],
      autod2 = Matrix::diag(autoDM)[Sgamma1]
    ))[,1:11]
# should be denseDiag
    (bob[1,]/2 + bob[2,])[1:5]


    thirdDiag[thirdDiag$i == Dorig,][1:5,]
    resThird$denseDiag[1:10,Dorig]


    thirdDiag[which.min(abs(thirdDiag$x - denseTkii[2,5])), ]


    bob = thirdList[[Dorig]]
    bob2 = bob[Sgamma1, Sgamma1]
    bob2[1:12,1:10]
    D3[1:12,1:10]

#  (D3/as.matrix(autoD))[1:5,1:5]
  } # end testing

#  Determinant
# to do: compute in data frame.  merge and sum


  Strace = unlist(lapply(thirdList, function(Hinv, Tuux, Sgamma) {
      sum(Matrix::diag(Hinv %*% Tuux[Sgamma,Sgamma]))
#    if(is.null(Tuux)) return(0)
#      sum((Hinv * Tuux[Sgamma,Sgamma])@x)  
  },  
  Hinv=invHessian, Sgamma = Sgamma1))
  Sparams = setdiff(seq(1, length(config$sparsity$full$p)), Sgamma1)

  StraceW = Strace[Sparams]
  StraceV = Strace[Sgamma1]
  HuuHutheta = invHessian %*% fullHessian[Sgamma1, Sparams]
  VHH = Matrix::drop(StraceV %*% HuuHutheta)

  result$extra = data.frame(
    param = Sparams,
    theta = resThird$first[-Sgamma1],
    W = StraceW,
    VHH = VHH,
    GHH =  Matrix::drop(resThird$first[Sgamma1] %*% HuuHutheta) 
  )
  result$extra$grad = result$extra$theta + result$extra$W + result$extra$VHH

  result$gradL = result$extra$grad

  if(all(deriv == 1)) {
    return(result$gradL)
  } 
  result$wrappers = wrappers_gamma
  result$config = config
  result$first = resThird$first


  result$third = resThird
  result$third$thirdList = thirdList

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
