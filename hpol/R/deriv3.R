#' @export
thirdDeriv = function(x, data, config, adFun, extra=FALSE) {


   # computing T_kii and H_ki, columnns are i, rows are k 
  # library("hpolcc");x=res$parameters_and_gamma;config=res$config;data=res$tmb_data;config$verbose=TRUE;config$dense=FALSE;config$num_threads = 10


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

  Nparameters = length(x)
  Spars0 = seq(0L, length=Nparameters)
  Sgamma0 = seq(Nbeta, len=Ngamma)
  Sgamma1 = Sgamma0+1L
  SbetaTheta0 = setdiff(Spars0, Sgamma0)

  if(missing(adFun)) {
    if(missing(data)) {
      stop("provide data if adFun is missing")
    }
    adFun = getAdFun(x, data=data, config=config)
  }

  resThird = thirdStrata(x, data, config, adFun) 
  
  if(identical(config$dense, TRUE)) {

    diagMat = matrix(resThird$diag, Nparameters)
    thirdDiag = data.frame(k=c(row(diagMat))-1, j=c(col(diagMat))-1, Tijk = c(diagMat))
    thirdDiag$i= thirdDiag$j

    fullHessian = matrix(resThird$second, Nparameters)

  } else { # sparse

  thirdDiag = data.frame(
    i=config$sparsity$second$nonSymmetric$j, 
    j=config$sparsity$second$nonSymmetric$j, 
    k = config$sparsity$second$nonSymmetric$i,
    Tijk = resThird$diag
  )

  fullHessian = Matrix::sparseMatrix(
    i = config$sparsity$second$full$i, 
    j = config$sparsity$second$full$j,
    x=resThird$second,
    dims = rep(length(x), 2), symmetric=TRUE, index1=FALSE)

}

cholHessianRandom = Matrix::Cholesky(fullHessian[Sgamma1, Sgamma1])
invHessianRandom = Matrix::solve(cholHessianRandom)
dUhat = as.matrix(- invHessianRandom %*% fullHessian[Sgamma1, -Sgamma1])

  thirdTensorAgg = as.data.frame(
    config$sparsity$third$ijk
  )
  thirdTensorAgg$Tijk = resThird$Tijk

thirdList = mapply(
  thirdTensorFromIndex,
  thirdIndexOffDiag = config$sparsity$third$index$offDiag, 
  thirdIndexDiag= config$sparsity$third$index$diag,
  MoreArgs = list(Tijk = thirdTensorAgg, Tiij=thirdDiag, N=Nparameters), 
  SIMPLIFY=FALSE 
)

dHlist = mapply(
  forDhList,
  Tijk = thirdList[Sgamma1],
  dUp = as.list(as.data.frame(t(as.matrix(dUhat)))),
  MoreArgs = list(Sgamma1=Sgamma1), SIMPLIFY=FALSE
)

dHlong <- data.table::rbindlist(dHlist, use.names = TRUE, fill = TRUE)
data.table::setnames(dHlong,
         old = names(dHlong)[-(1:2)],
         new = paste0("p", SbetaTheta0))

# fast aggregate: sum all columns except i,j by i,j
dHagg <- dHlong[, lapply(.SD, sum, na.rm = TRUE), by = .(i, j)]

dH = mapply(forDh, 
  Dpar = SbetaTheta0,
  MoreArgs = list(dHagg = dHagg, thirdList = thirdList,
    Sgamma1 = Sgamma1, dims=c(Ngamma, Ngamma)),
  SIMPLIFY=FALSE 
)

dDet = mapply(traceProd, Dhp = dH, 
  MoreArgs = list(Hinv = invHessianRandom), 
  SIMPLIFY=TRUE
) / 2


result = list(  
  halfLogDet = drop(Matrix::determinant(cholHessianRandom, log=TRUE, sqrt=TRUE)$modulus),
  dUhat = dUhat, dDet = dDet,    dH = dH, 
  fullHessian=fullHessian, invHessianRandom = invHessianRandom,  first = resThird$first
)

if(identical(extra, TRUE)) {
  result = c(result, list(
    thirdList = thirdList, raw = resThird))
}

return(result)
}
