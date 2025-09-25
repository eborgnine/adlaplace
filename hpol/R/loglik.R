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
  if(missing(gamma_start)) {
    gamma_start = rep(0, Ngamma)
  }

#library('hpolcc')
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

   # computing T_kii and H_ki, columnns are i, rows are k 
#library('hpolcc')
    configDiag = config;configDiag$dense = TRUE
    resThirdDiag = thirdDiagonals(
      result$fullParameters, data, configDiag
    ) 
    resThirdOffDiag = thirdOffDiagonals(
      result$fullParameters, data, config
    ) 

    if(identical(configDiag$dense, TRUE)) {
          # T_iik, doubles are columns of resThirdDiag$diag
      fullHessian = Matrix::Matrix(resThirdDiag$second)
      thirdDiag = data.frame(
        k = rep(seq(0, len=ncol(resThirdDiag$diag)), each=nrow(resThirdDiag$diag)),
        i = rep(seq(0, len=nrow(resThirdDiag$diag)), ncol(resThirdDiag$diag)),
        x = 2*as.vector(resThirdDiag$diag)
      )
      thirdDiag$j = thirdDiag$i
    } else {
    thirdDiag = data.frame(
      i = config$sparsity$second$nonSymmetric$i,
      j = config$sparsity$second$nonSymmetric$i,
      k = config$sparsity$second$nonSymmetric$j,
      x = 2*(resThirdDiag$diag)
    )
        fullHessian = Matrix::forceSymmetric(
      Matrix::sparseMatrix(
      j = config$sparsity$second$nonSymmetric$j,
      i = config$sparsity$second$nonSymmetric$i,
      x = drop(resThirdDiag$second),
      dims = rep(length(resThirdDiag$first), 2), index1=FALSE
    ))

    }

    if(identical(config$dense, TRUE)) {

      thirdNonDiag = config$sparsity$third$pairs[
        rep(1:nrow(config$sparsity$third$pairs), nrow(resThirdOffDiag)),c('i','j')]
      thirdNonDiag$k = rep(seq(0, len=nrow(resThirdOffDiag)), ncol(resThirdOffDiag))
      thirdNonDiag$taylor3 = as.vector(resThirdOffDiag)
    } else {
      thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
      thirdNonDiag$taylor3 = drop(resThirdOffDiag)
    }

  # 3rd taylor is   T_iik/2 + T_jjk/2 + T_ijk
  # pairs are i and j
    thirdDiag$ik = apply(thirdDiag[,c('i','k')],1,paste,collapse='_')
    thirdNonDiag$ik = apply(thirdNonDiag[,c('i','k')],1,paste,collapse='_')
    thirdNonDiag$jk = apply(thirdNonDiag[,c('j','k')],1,paste,collapse='_')
    matchIik = match(thirdNonDiag$ik, thirdDiag$ik)
    matchJjk = match(thirdNonDiag$jk, thirdDiag$ik)
    thirdNonDiag$Tiik = thirdDiag[matchIik, 'x']
    thirdNonDiag$Tjjk = thirdDiag[matchJjk, 'x']

    thirdNonDiag$x = thirdNonDiag$taylor3 - 0.375*(
      thirdNonDiag$Tiik  + thirdNonDiag$Tjjk)


    theCols = c('i','j','k', 'x')
    third = rbind(
      thirdDiag[,theCols],
      thirdNonDiag[,theCols]
    )

    thirdList = mapply(
      function(k, third, N) {
        thirdHere = third[apply(third[,c('i','j','k')] == k, 1, any), ]
        if(!nrow(thirdHere)) return(NULL)
          newxy = t(apply(thirdHere[,c('i','j','k')], 1, function(xx) sort(c(xx[xx!=k], rep(k,2))[1:2] ) ))
        try(Matrix::sparseMatrix(i=newxy[,1], j=newxy[,2], x=thirdHere[,'x'],
          dims = rep(N, 2), symmetric=TRUE, index1=FALSE))
      },
      k = seq(from=0, len=nrow(fullHessian)), 
      MoreArgs = list(third=third, N=nrow(fullHessian))
    )

  cholHessianRandom = Matrix::Cholesky(fullHessian[Sgamma1, Sgamma1])
  invHessianRandom = Matrix::solve(cholHessianRandom)
  result$logDetHessian = drop(Matrix::determinant(
      cholHessianRandom, log=TRUE, sqrt=FALSE
    )$modulus)

    result$minusLogLik = 
      result$fval +
      as.numeric(result$logDetHessian)/2 + 
      0.5 * Ngamma * 1.8378770664093454835606594728  # log 2 pi


  Strace = unlist(lapply(thirdList, function(Hinv, Tuux, Sgamma1) {
      sum(Matrix::diag(Hinv %*% Tuux[Sgamma1,Sgamma1]))
#    if(is.null(Tuux)) return(0)
#      sum(as(Hinv * Tuux[Sgamma1,Sgamma1], "generalMatrix")@x)  
  },  
  Hinv=invHessianRandom, Sgamma1 = Sgamma1))

  Sparams = setdiff(1:nrow(fullHessian), Sgamma1)

  StraceW = Strace[Sparams]
  StraceV = Strace[Sgamma1]
  HuuHutheta = invHessianRandom %*% fullHessian[Sgamma1, Sparams]
  VHH = Matrix::drop(StraceV %*% HuuHutheta)

  result$extra = data.frame(
    param = Sparams,
    theta = resThirdDiag$first[-Sgamma1,1],
    W = StraceW,
    VHH = VHH,
    GHH =  Matrix::drop(resThirdDiag$first[Sgamma1] %*% HuuHutheta) 
  )

  result$extra$grad = result$extra$theta + result$extra$W + result$extra$VHH
  result$gradL = result$extra$grad

  if(all(deriv == 1)) {
    return(result$gradL)
  } 

  result$wrappers = wrappers_gamma
  result$config = config
  result$first = resThirdDiag$first[,1]
  result$fullHessian = fullHessian

    result$third = thirdList
    result$diag = resThirdDiag
    result$nonDiag = thirdNonDiag

  return(result)
}
