#' @export

sparsity_grouped = function(x, data, config, verbose=FALSE) {

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
	dataNoMap$Qdiag = rep(1, 0)

	if(verbose) cat("getting first deriv...")
	firstDeriv = gradLogical(x, dataNoMap, config)
	if(verbose) cat("done\ngetting clusters...")
	km <- try(kmeans(t(firstDeriv), centers = ceiling(1.1*Nclusters), 
		iter.max=25000, nstart=10*Nclusters, algorithm='Hartigan-Wong'))
	if(verbose) cat("done\n")

	if(!any(class('km') == 'try-error')) {
		mergeThreshold = floor(0.25*Nstrata/Nclusters)
		NtoMerge = (which(cumsum(sort(km$size)) < mergeThreshold))
		if(length(NtoMerge)) {
			NtoMerge = max(NtoMerge)
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
		}
	} else {
		pr = prcomp(firstDeriv)
		theOrder = order(pr$rotation[,1])
		theCl = floor(seq(1, Nclusters+0.999, len= Nstrata))
		km = list(cluster = theCl)
	}


	clusterTable = table(km$cluster)
	tableOrder = order(clusterTable, decreasing=TRUE)
	km$cluster = match(km$cluster, tableOrder)


	firstDeriv = Matrix::Matrix(firstDeriv)

	strataMatrix = Matrix::sparseMatrix(
		i = seq_len(Nstrata),
		j = km$cluster,
		dims = c(Nstrata, max(km$cluster))
	)

	strataMatrixList = list(i=strataMatrix@i, j = as(strataMatrix, 'TsparseMatrix')@j, p=strataMatrix@p)

	strataListForHessian = mapply(
		function(x, y) {
			list(list(p = y$p[c(x, x+1)], i=y$i))
		},
		x = seq(1,len=ncol(strataMatrix)), 
		MoreArgs = list(y=strataMatrixList))

#	hessianDenseAll = hessianDenseLogical(parameters=x, data=dataNoMap,config=config, strata = res$groups$groups)
	# now find hessian for each block
	if(verbose) cat("getting hessian by group...")
	hessianByBlock = parallel::mcmapply(function(strata, config, ...) {
		hessianDenseLogical(..., config = c(config, list(groups = strata)))
	},
		strata = strataListForHessian,
		MoreArgs = list(
			parameters=x, data=dataNoMap, 
			config = c(config[setdiff(names(config), 'num_threads')], 
				list(num_threads=1))),
		SIMPLIFY=FALSE, mc.cores=config$num_threads
	)
	hessianByBlock2 = lapply(hessianByBlock, Matrix::Matrix, sparse=TRUE)

	#hessian for random part,
	if(verbose) cat("getting Q hessian")
	hessianQ = list(dense=hessianQdense(parameters=x, data=data, config=config))


# get non-zeros of tensor

	hessianQ$hessian = Matrix::forceSymmetric(Matrix::Matrix(hessianQ$dense, sparse=TRUE))
	hessianQT = as(as(hessianQ$hessian, 'generalMatrix'),'TsparseMatrix')
	hessianQns = as(hessianQ$hessian, 'generalMatrix')
	if(verbose) cat("getting third sparsity...")
	ijkQ = getThirdFromHessian(hessianQ$hessian)
	if(verbose) cat("done\n")
	if(any(ijkQ[,'Nunique'] == 3)) warning("need to implement non-diagonal Q third")

	fullHessian=do.call(rbind, lapply(hessianByBlock2, function(xx) cbind(xx@i,as(xx, "TsparseMatrix")@j)))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = t(apply(fullHessian, 1, sort))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = fullHessian[order(fullHessian[,2], fullHessian[,1]),]

	hessianNoQ = Matrix::sparseMatrix(i=apply(fullHessian,1,min), j=apply(fullHessian,1,max), 
		x=rep(1, nrow(fullHessian)), symmetric=TRUE, index1=FALSE)

	fullHessian = rbind(fullHessian, cbind(hessianQT@i, hessianQT@j))
	fullHessian = t(apply(fullHessian, 1, sort))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = fullHessian[order(fullHessian[,2], fullHessian[,1]),]

	fullHessianMatrix = Matrix::sparseMatrix(i=fullHessian[,1], j=fullHessian[,2], symmetric=TRUE, index1=FALSE,
		dims = c(Ntotal, Ntotal), repr='T')
	if(verbose) cat("getting optimal pairs")
	fullList = getOptimalPairs(fullHessianMatrix, Sparams=Sparams, Sgamma1=Sgamma1)
	if(verbose) cat("done\n")


	fullHessianPairs = paste(fullList$second$full$i, fullList$second$full$j, sep='_')
	fullHessianPairsNs = paste(fullList$second$nonSymmetric$i,fullList$second$nonSymmetric$j, sep='_')
	fullHessianPairsR = paste(fullList$second$random$i, fullList$second$random$j, sep='_')
	fullHessianPairsRNs = paste(fullList$second$random$i, fullList$second$random$j, sep='_')

	hessianQrandom = hessianQ$hessian[Sgamma1, Sgamma1] 
	hessianQrandomNs = as(hessianQrandom, 'generalMatrix')

	sparsityQ = list(
		full=list(i=hessianQ$hessian@i, p=hessianQ$hessian@p, j = as(hessianQ$hessian, 'TsparseMatrix')@j),
		nonSymmetric=list(i=hessianQns@i, p=hessianQns@p, j = as(hessianQns, 'TsparseMatrix')@j),
		random = list(i=hessianQrandom@i, p=hessianQrandom@p, j = as(hessianQrandom, 'TsparseMatrix')@j),
		randomNS = list(i=hessianQrandomNs@i, p=hessianQrandomNs@p, j = as(hessianQrandomNs, 'TsparseMatrix')@j)
	)

	sparsityQ$full$match = try(match(
		paste(sparsityQ$full$i, sparsityQ$full$j, sep='_'), 
		fullHessianPairs
	)-1L)
	sparsityQ$random$match = try(match(
		paste(sparsityQ$random$i, sparsityQ$random$j, sep='_'), 
		fullHessianPairsR
	)-1L)


	# find full hessian sparsity
	# for each strata, get index in full hessian
	if(verbose) cat("getting sparsity by block...")
	sparsityList = parallel::mcmapply(getOptimalPairs,
		hessian = hessianByBlock2,
		MoreArgs = list(Sparams = Sparams, Sgamma1=Sgamma1, 
			hessianPairs = fullHessianPairs,
			hessianPairsR = fullHessianPairsR), 
		SIMPLIFY=FALSE, mc.cores=config$num_threads)
	if(verbose) cat("done\n")

	
	result = list(
		group_sparsity = sparsityList,
		groups = strataMatrixList,
		sparsity = c(fullList, list(Q=sparsityQ)),
		firstDeriv = firstDeriv
	)

}