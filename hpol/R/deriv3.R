#' @export
thirdDeriv = function(x, data, config) {

  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Nparameters = length(x)
  Spars0 = seq(0, length=Nparameters)
  Sgamma0 = seq(Nbeta, len=Ngamma)
  Sgamma1 = Sgamma0+1
  SbetaTheta0 = setdiff(Spars0, Sgamma0)

   # computing T_kii and H_ki, columnns are i, rows are k 
  resThirdDiag = thirdDiagonals(
    x, data, config
  ) 

  resThirdOffDiag = thirdOffDiagonals(
    x, data, config
  ) 

  if(identical(config$dense, TRUE)) {
          # T_iik, doubles are columns of resThirdDiag$diag
    fullHessian = Matrix::Matrix(resThirdDiag$second)
    thirdDiag = data.frame(
      i = rep(seq(0, len=nrow(resThirdDiag$diag)), each=ncol(resThirdDiag$diag)),
      k = rep(seq(0, len=ncol(resThirdDiag$diag)), nrow(resThirdDiag$diag)),
      x = 2*as.vector(resThirdDiag$diag)
    )

    thirdNonDiag = config$sparsity$third$pairs[
    	rep(1:nrow(config$sparsity$third$pairs), each=nrow(resThirdOffDiag)),
    	c('i','j')]
    thirdNonDiag$k = rep(seq(0, len=nrow(resThirdOffDiag)), ncol(resThirdOffDiag))
    thirdNonDiag$taylor3 = as.vector(resThirdOffDiag)
    thirdNonDiag = thirdNonDiag[
    	apply(
    		thirdNonDiag[,c('i','j','k')], 1, lengthUnique
    	)==3, ]
  } else { # sparse
    fullHessian = Matrix::forceSymmetric(
      Matrix::sparseMatrix(
        i = config$sparsity$second$nonSymmetric$i,
        j = config$sparsity$second$nonSymmetric$j,
        x = drop(resThirdDiag$second),
        dims = rep(length(resThirdDiag$first), 2), index1=FALSE
      ))

    thirdDiag = data.frame(
      i = config$sparsity$second$nonSymmetric$j,
      k = config$sparsity$second$nonSymmetric$i,
      x = 2*as.vector(resThirdDiag$diag)
    )

    thirdNonDiag = config$sparsity$third$ijk[,c('i','j','k')]
    thirdNonDiag$taylor3 = drop(resThirdOffDiag)
  }
    thirdDiag$j = thirdDiag$i


  # 3rd taylor is   T_iik/2 + T_jjk/2 + T_ijk
  # pairs are i and j
  thirdDiag$ik = apply(thirdDiag[,c('i','k')],1,paste,collapse='_')
  thirdNonDiag$ik = apply(thirdNonDiag[,c('i','k')],1,paste,collapse='_')
  thirdNonDiag$jk = apply(thirdNonDiag[,c('j','k')],1,paste,collapse='_')
  matchIik = match(thirdNonDiag$ik, thirdDiag$ik)
  matchJjk = match(thirdNonDiag$jk, thirdDiag$ik)
  thirdNonDiag$Tiik = thirdDiag[matchIik, 'x']
  thirdNonDiag$Tjjk = thirdDiag[matchJjk, 'x']

  thirdNonDiag$x = thirdNonDiag$taylor3 - 0.5*(thirdNonDiag$Tiik  + thirdNonDiag$Tjjk)
  
# might be some duplicated
  nonDiagUnique = t(apply(thirdNonDiag[,c('i','j','k')], 1, sort))  
  thirdNonDiag = thirdNonDiag[!duplicated(nonDiagUnique), ]


  theCols = c('i','j','k', 'x')
  third = rbind(
    thirdDiag[abs(thirdDiag$x) > 1e-20,theCols],
    thirdNonDiag[abs(thirdNonDiag$x) > 1e-20,theCols]
  )
  third = third[order(third[,'i'], third[,'j'], third[,'k']),]

  # to do: don't cmpute third list, cmpute dHlist without it.

   
  cholHessianRandom = Matrix::Cholesky(fullHessian[Sgamma1, Sgamma1])
  invHessianRandom = Matrix::solve(cholHessianRandom)
  dUhat = - invHessianRandom %*% fullHessian[Sgamma1, -Sgamma1] 


  thirdList = mapply(
  	thirdTensor,
  	k=seq(from=0, len=Nparameters),
  	MoreArgs = list(third=third, N=Nparameters)
  )


  dHlist = mapply(
  	function(Tijk, dUp, Sgamma1) {
      if(is.null(Tijk)) return(NULL)
  		There =try( as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix'))
  		cbind(j=There@i, i=There@j, outer(There@x, dUp))
  	},
  	Tijk = thirdList[Sgamma1],
  	dUp = as.list(as.data.frame(t(as.matrix(dUhat)))),
  	MoreArgs = list(Sgamma1=Sgamma1)
  )

  TijpAdd = mapply(function(Tijk, Sgamma1) {
  if(is.null(Tijk)) return(NULL)
	There = as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix')
	cbind(i=There@i, j=There@j, x=There@x)
  }, Tijk = thirdList[-Sgamma1], MoreArgs = list(Sgamma1=Sgamma1))
  names(TijpAdd) = paste0("p", SbetaTheta0)

# try to do without thirdList
  if(FALSE) {
  dHlist = mapply(
  	function(Dgamma, dUp, third, Sgamma1,  Nparameters) {
  		Tijk = thirdTensor(k=Dgamma, third=third, N=Nparameters)
  		There = as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix')
  		cbind(j=There@i, i=There@j, outer(There@x, dUp))
  	},
  	Dgamma = Sgamma0,
  	dUp = as.list(as.data.frame(t(as.matrix(dUhat)))),
  	MoreArgs = list(Sgamma1=Sgamma1, third = third, Nparameters=Nparameters)
  )
  TijpAdd = mapply(function(Dpar, third, Sgamma1,  Nparameters) {
  	Tijk = thirdTensor(k=Dpar, third=third, N=Nparameters)  	
	There = as(Tijk[Sgamma1,Sgamma1], 'TsparseMatrix')
	cbind(i=There@i, j=There@j, x=There@x)
  }, Dpar = SbetaTheta0, MoreArgs = list(third=third, Sgamma1=Sgamma1, Nparameters=Nparameters))
  names(TijpAdd) = paste0("p", SbetaTheta0)

}


  dHlong = as.data.frame(do.call(rbind, dHlist))
  colnames(dHlong) = c(names(dHlong)[1:2], paste0("p", SbetaTheta0))
  dHagg = aggregate(
  	dHlong[,setdiff(names(dHlong), c('i','j')), drop=FALSE], 
  	dHlong[,c('i','j')], 
  	sum)

#  dHlong$k = rep(1:length(dHlist), unlist(lapply(dHlist, nrow)))


  
  dH = mapply(function(TU, ij, Tijp, dims) {
  	toAgg = rbind(cbind(ij, x=TU), Tijp)
  	theAgg = aggregate(toAgg[,'x', drop=FALSE], toAgg[,c('i','j')], sum, na.rm=TRUE)
  	Matrix::sparseMatrix(
  		i=pmax(theAgg$j, theAgg$i), 
  		j=pmin(theAgg$i, theAgg$j), x=theAgg$x, 
  		index1=FALSE, dims=dims, symmetric=TRUE)
  }, 
  TU = as.list(dHagg[,names(TijpAdd)]),
  Tijp = TijpAdd,
  MoreArgs = list(ij = dHagg[,c('i','j')], dims=c(Ngamma, Ngamma))
	)


  # i,j are gamma indices


dDet = mapply(function(Dhp, Hinv) {
	sum(Matrix::diag(Hinv %*% Dhp))
	sum((Hinv * Dhp))
}, Dhp = dH, MoreArgs = list(Hinv = invHessianRandom)
)


  if(FALSE) {
  	 pair = pairsGamma[1,]# go to getDh
  	Dp = 2
  	DHpInitial = as.matrix(thirdList[[Dp]][Sgamma1, Sgamma1])
  	DHp1 = 0
  	bob=rep(NA, length(Sgamma1))
  	for(D in 1:length(Sgamma1)) {
  		DHp1 = DHp1 + as.matrix(thirdList[[Sgamma1[D] ]][Sgamma1,Sgamma1])*dUhat[D,Dp]
  	}
#  	stuff1 = unlist(lapply(thirdList[Sgamma1], function(xx) xx[DparL2, DparL3]))
sum(dUhat[,2] * as.vector(thirdHere[Sgamma1]))
DHp1[1:25,1:7]
sum(bob*dUhat[,Dp])
  	DHp = DHp1 + DHpInitial
DHp[1:15,1:5]
dHlist[[Dp]][1:15,1:5]

dH[[2]][1:3,1:3]
DHp[1:3,1:3]
stuff = DHp - dH[[2]]

  }



list(fullHessian=fullHessian, third =third, thirdList = thirdList,  first = resThirdDiag$first, 
	invHessianRandom = invHessianRandom,
	halfLogDet = drop(Matrix::determinant(cholHessianRandom, log=TRUE, sqrt=TRUE)$modulus),
	dUhat = dUhat, dH = dH, dDet = dDet
)

}