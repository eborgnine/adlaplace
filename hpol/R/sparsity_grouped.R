#' @export
 
sparsity_grouped = function(x, data, config) {

# library('hpolcc');x = res$parameters_for_sparsity;data=res$tmb_data;config=res$config

	Nstrata = ncol(data$cc_matrixTp)
	Sstrata = seq(0, len=Nstrata)
	Nclusters  <- config$num_threads*2                                   

  Nbeta = nrow(data$XTp)
  Ngamma = nrow(data$ATp)
  Ntotal = length(x)
  Sgamma1 = seq(from=Nbeta+1, len=Ngamma)
  Nparams = Ntotal - Ngamma  
  Sparams = setdiff(seq(0, len=Ntotal), seq(from=Nbeta, len=Ngamma))


	x[x==0] = 0.01

	singleStrataList = mapply(function(x) {
		list(p=c(0,1), i=x)
	}, x=seq(0,len=Nstrata), SIMPLIFY=FALSE)

	# get rid of map so Q won't be added to likelihood
	dataNoMap = data


	firstDeriv = mapply(grad,
		strata = singleStrataList,
		MoreArgs = list(
			parameters=x, data=dataNoMap, sparsity = list(),
			config = c(config[
				setdiff(names(config), c('num_threads'))
				], list(num_threads=1))),
		SIMPLIFY=FALSE
	)
	firstDeriv = do.call(cbind, lapply(firstDeriv, function(xx) {xx$grad > 0} ))

		km <- try(kmeans(t(firstDeriv), centers = ceiling(1.5*Nclusters), iter.max=15000,nstart=5*Nclusters, algorithm='Hartigan-Wong'))
		if(!any(class('km') == 'try-error')) {
			mergeThreshold = floor(0.5*Nstrata/Nclusters)
			NtoMerge = max(which(cumsum(sort(km$size)) < mergeThreshold))
			mergeSize = sort(km$size)[NtoMerge]
			whichToMerge = which(km$size <= mergeSize)

			km$cluster2 = km$cluster
			km$cluster2[km$cluster %in% whichToMerge] = min(whichToMerge)
			km$cluster = as.integer(factor(km$cluster2))

			# split biggest cluster
			theTable = table(km$cluster)
			if(max(theTable) > 1.5*order(theTable, decreasing=TRUE)[2]) {
				whichInBiggest = which(km$cluster == which.max(theTable))
				firstDeriv2 = firstDeriv[,whichInBiggest]
				km2 = kmeans(t(firstDeriv2), centers = 2, iter.max=15000,nstart=200, algorithm='Hartigan-Wong')
				km$cluster[whichInBiggest] = -km2$cluster
				km$cluster = as.integer(factor(km$cluster))
			}

		} else {
			pr = prcomp(firstDeriv)
			theOrder = order(pr$rotation[,1])
			theCl = floor(seq(1, Nclusters+0.999, len= Nstrata))
			km = list(cluster = theCl)
	}

if(FALSE) {
	image(seq(0, len=nrow(firstDeriv)), seq(0, len=ncol(firstDeriv)), firstDeriv[,order(km$cluster)])
	abline(h=cumsum(table(km$cluster)), col='blue', lty=1)

	length(unique(km$cluster))
	sort(table(km$cluster))
}

	strataMatrix = Matrix::sparseMatrix(
		i = seq_len(Nstrata),
		j = km$cluster,
		dims = c(Nstrata, max(km$cluster))
	)
	strataMatrix = strataMatrix[,order(table(km$cluster), decreasing=TRUE)]

	strataMatrixList = list(i=strataMatrix@i, j = as(strataMatrix, 'TsparseMatrix')@j, p=strataMatrix@p)

	strataListForHessian = mapply(
		function(x, y) {
			list(list(p = y$p[c(x, x+1)], i=y$i))
		},
		x = seq(1,len=ncol(strataMatrix)), 
		MoreArgs = list(y=strataMatrixList))


	# now find hessian for each block
	hessianByBlock = parallel::mcmapply(hessianDense,
		strata = strataListForHessian,
		MoreArgs = list(
		parameters=x, data=dataNoMap, sparsity = list(),
		config = c(config[setdiff(names(config), 'num_threads')], 
			list(num_threads=1))),
		SIMPLIFY=FALSE, mc.cores=config$num_threads
	)
	hessianByBlock2 = lapply(hessianByBlock, function(xx) Matrix::Matrix(xx$hessian, sparse=TRUE))

	#hessian for random part,
	hessianQ = list(dense=hessianQdense(parameters=x, data=data, config=config))


# get non-zeros of tensor

	hessianQ$hessian = Matrix::forceSymmetric(Matrix::Matrix(hessianQ$dense, sparse=TRUE))
	hessianQT = as(as(hessianQ$hessian, 'generalMatrix'),'TsparseMatrix')
	hessianQns = as(hessianQ$hessian, 'generalMatrix')
    ijkQ = getThirdFromHessian(hessianQ$hessian)
    if(any(ijkQ[,'Nunique'] == 3)) warning("need to implement non-diagonal Q third")

	fullHessian=do.call(rbind, lapply(hessianByBlock2, function(xx) cbind(xx@i,as(xx, "TsparseMatrix")@j)))
	fullHessian = fullHessian[!duplicated(fullHessian), ]

	fullHessian = rbind(fullHessian, cbind(hessianQT@i, hessianQT@j))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = fullHessian[fullHessian[,1] <= fullHessian[,2], ]
	fullHessianMatrix = Matrix::sparseMatrix(i=fullHessian[,1], j=fullHessian[,2], symmetric=TRUE, index1=FALSE,
		dims = c(Ntotal, Ntotal), repr='T')
	fullList = getOptimalPairs(fullHessianMatrix, Sparams=Sparams, Sgamma1=Sgamma1)


  	fullHessianPairs = paste(fullList$second$full$i, fullList$second$full$j, sep='_')
  	fullHessianPairs2 = paste(fullList$second$nonSymmetric$i,fullList$second$nonSymmetric$j, sep='_')
  	fullHessianPairsR = paste(fullList$second$random$i, fullList$second$random$j, sep='_')

  	randomFromFull = paste(fullList$second$random$i+Nbeta, fullList$second$random$j+Nbeta, sep='_')


	sparsityQ = list(
		full=list(i=hessianQ$hessian@i, p=hessianQ$hessian@p, j = as(hessianQ$hessian, 'TsparseMatrix')@j),
		nonSymmetric=list(i=hessianQns@i, p=hessianQns@p, j = as(hessianQns, 'TsparseMatrix')@j)
	)

  	sparsityQ$full$match = try(match(
  		paste(sparsityQ$full$i, sparsityQ$full$j, sep='_'), 
  		fullHessianPairs
  	))
  	sparsityQ$nonSymmetric$match = try(match(
  		paste(sparsityQ$nonSymmetric$i, sparsityQ$nonSymmetric$j, sep='_'), 
  		fullHessianPairs2
  	))


	# find full hessian sparsity
	# for each strata, get index in full hessian

	sparsityList = parallel::mcmapply(getOptimalPairs,
		hessian = hessianByBlock2,
		MoreArgs = list(Sparams = Sparams, Sgamma1=Sgamma1, 
			hessianPairs = fullHessianPairs,
			hessianPairsNS = fullHessianPairs2,
			hessianPairsR = fullHessianPairsR,
			randomFromFull = randomFromFull), 
		SIMPLIFY=FALSE, mc.cores=config$num_threads)

	
	result = list(
		group_sparsity = sparsityList,
		groups = strataMatrixList,
		sparsity = c(fullList, list(Q=sparsityQ))
	)

}