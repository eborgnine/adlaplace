#' @export
thirdDeriv = function(x, data, config, adFun, extra=FALSE) {


   # computing T_kii and H_ki, columnns are i, rows are k 
  # library("hpolcc");x=res$parameters_and_gamma;config=res$config;data=res$tmb_data;config$verbose=FALSE;config$dense=TRUE#;config$num_threads = 10
  if(missing(adFun)) {
    if(missing(data)) stop("provide data if adFun is missing")
    adFun = getadFun(x, data=data, config=config)
  }

  Nparameters = length(x)
  Spars0 = seq(0, length=Nparameters)

  if(!missing(data)) {
    Nbeta = nrow(data$XTp)
    Ngamma = nrow(data$ATp)
  } else {
    Nbeta = length(config$beta)
    Ntheta = length(config$theta)
    Ngamma = Nparameters - Nbeta - Ntheta
    if(Nbeta == 0) warning("assuming no betas.  add data or config$beta if this is incorrect")
    if(Ntheta == 0) warning("assuming no betas.  add data or config$beta if this is incorrect")
  }

  Sgamma0 = seq(Nbeta, len=Ngamma)
  Sgamma1 = Sgamma0+1
  SbetaTheta0 = setdiff(Spars0, Sgamma0)


  resThird = thirdStrata(
    x, data, config, adFun
  ) 

  
  if(identical(config$dense, TRUE)) {


    thirdTensorList = mapply(function(group_sparsity, Tijk)
      {
        Npairs = nrow(group_sparsity$third$pairs)
        resultIJK = as.data.frame(t(apply(cbind(group_sparsity$third$pairs[
            rep(1:Npairs, each=Nparameters), c('i','j')],
          k=rep(seq(0, len=Nparameters), Npairs)), 1, sort)))
          colnames(resultIJK) = c('i','j','k')
          duplicatedHere = duplicated(resultIJK)
          resultIJK$Tijk = Tijk
        result = resultIJK[!duplicatedHere,]
        result = result[result$Tijk != 0, ]        
      }, group_sparsity = config$group_sparsity, Tijk = resThird$Tijk,
      SIMPLIFY=FALSE
    )

  thirdDiag = expand.grid(i=seq(0, len=Nparameters), k=seq(0, len=Nparameters))
  thirdDiag$j = thirdDiag$i
  thirdDiag$Tijk = as.vector(resThird$diag)

  fullHessian = matrix(resThird$second, Nparameters)

  } else {

  thirdTensorList = mapply(function(third, sparsity) {
      ijkHere = sparsity$third$ijk[,c('i','j','k')]
      ijkHere$Tijk = third
      ijkHere
  },
  third = resThird$Tijk[1:10],
  sparsity = config$group_sparsity[1:10],
  SIMPLIFY=FALSE)

  thirdDiag = data.frame(
    i=config$sparsity$second$nonSymmetric$i, 
    j=config$sparsity$second$nonSymmetric$i, 
    k = config$sparsity$second$nonSymmetric$j,
    Tijk = resThird$diag
  )

  fullHessian = Matrix::sparseMatrix(
    i = config$sparsity$second$full$i, 
    j = config$sparsity$second$full$j,
    x=resThird$second,
    dims = rep(Nparameters, 2), symmetric=FALSE, index1=FALSE)

  }


  thirdTensorDf = do.call(rbind, thirdTensorList)
  thirdTensorDf = cbind(as.data.frame(t(apply(thirdTensorDf[,c('i','j','k')], 1, sort))), thirdTensorDf[,'Tijk',drop=FALSE])
  colnames(thirdTensorDf) = c('i','j','k','Tijk')

  thirdTensorAgg = aggregate(
    thirdTensorDf[,'Tijk', drop=FALSE], 
    as.data.frame(thirdTensorDf[,c('i','j','k')]), 
    sum)

  # get rid of diagonals, only has an effect if dense=TRUE
  Nunique = apply(thirdTensorAgg[,c('i','j','k')], 1, lengthUnique)
  thirdTensorAgg = thirdTensorAgg[Nunique == 3, ]

  thirdAll = rbind(thirdTensorAgg, thirdDiag)
  thirdAll = thirdAll[order(thirdAll$i, thirdAll$j, thirdAll$k),]
  thirdAll = thirdAll[abs(thirdAll$Tijk)> .Machine$double.eps, ]
  colnames(thirdAll) = gsub("Tijk", "x", colnames(thirdAll))

  thirdList = parallel::mcmapply(
    thirdTensor,
    k=seq(from=0, len=Nparameters),
    MoreArgs = list(third=thirdAll, N=Nparameters), 
    mc.cores = pmax(config$num_threads, 1, na.rm=TRUE)
  )


  cholHessianRandom = Matrix::Cholesky(fullHessian[Sgamma1, Sgamma1])
  invHessianRandom = Matrix::solve(cholHessianRandom)
  dUhat = - invHessianRandom %*% fullHessian[Sgamma1, -Sgamma1] 

  dHlist = parallel::mcmapply(
    forDhList,
    Tijk = thirdList[Sgamma1],
    dUp = as.list(as.data.frame(t(as.matrix(dUhat)))),
    MoreArgs = list(Sgamma1=Sgamma1), 
    mc.cores = pmax(config$num_threads, 1, na.rm=TRUE)
  )

  TijpAdd = parallel::mcmapply(forTijpAdd, 
    Tijk = thirdList[-Sgamma1], MoreArgs = list(Sgamma1=Sgamma1),
    mc.cores = pmax(config$num_threads, 1, na.rm=TRUE)
  )
  names(TijpAdd) = paste0("p", SbetaTheta0)


  dHlong = as.data.frame(do.call(rbind, dHlist))
  colnames(dHlong) = c(names(dHlong)[1:2], paste0("p", SbetaTheta0))
  dHagg = aggregate(
    dHlong[,setdiff(names(dHlong), c('i','j')), drop=FALSE], 
    dHlong[,c('i','j')], 
    sum)

  
  dH = parallel::mcmapply(forDh, 
    TU = as.list(dHagg[,names(TijpAdd)]),
    Tijp = TijpAdd,
    MoreArgs = list(ij = dHagg[,c('i','j')], dims=c(Ngamma, Ngamma)), 
    mc.cores = pmax(config$num_threads, 1, na.rm=TRUE)
  )


dDet = parallel::mcmapply(traceProd, Dhp = dH, 
  MoreArgs = list(Hinv = invHessianRandom), 
    mc.cores = pmax(config$num_threads, 1, na.rm=TRUE)
  ) / 2


result = list(  
  halfLogDet = drop(Matrix::determinant(cholHessianRandom, log=TRUE, sqrt=TRUE)$modulus),
  dUhat = dUhat, dH = dH, dDet = dDet,
  fullHessian=fullHessian, invHessianRandom = invHessianRandom,  first = resThird$first
)

if(identical(extra, TRUE)) {
result = c(result, list(third =thirdAll, 
  thirdList = thirdList, raw = resThird))
}

return(result)
}
