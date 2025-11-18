#' @export

sparsity_grouped = function(x, data, config, verbose=FALSE) {

# library('hpolcc');x = res$parameters_and_gamma;data=res$tmb_data;config=res$config;verbose=TRUE


	Nstrata = ncol(data$cc_matrixTp)
	Sstrata = seq(0, len=Nstrata)

	Nclusters = config$Nclusters
	if(!length(Nclusters)) {
		Nclusters  <- config$num_threads*2                                   
	}

	Nbeta = nrow(data$XTp)
	Ngamma = nrow(data$ATp)
	Ntotal = length(x)
	Sgamma1 = seq(from=Nbeta+1, len=Ngamma)
	Sgamma0 = Sgamma1 - 1
	Nparams = Ntotal - Ngamma  
	Sparams = setdiff(seq(0, len=Ntotal), seq(from=Nbeta, len=Ngamma))


	singleStrataList = mapply(function(x) {
		list(p=c(0,1), i=x)
	}, x=seq(0,len=Nstrata), SIMPLIFY=FALSE)

	# get rid of map so Q won't be added to likelihood
	dataNoMap = data
	dataNoMap$Qdiag = rep(1, 0)

	if(verbose) {
		cat("getting first deriv...")
	}
	firstDeriv = gradLogical(x, dataNoMap, config)
	firstDeriv = Matrix::Matrix(firstDeriv)

	if(verbose) {
		cat("done\ngetting clusters...")
	}
	if(requireNamespace("RSpectra", quietly=TRUE) ) {
		e <- RSpectra::svds(1.1*firstDeriv, k = 1)
		loadings <- e$v		
	} else {
		X <- scale(firstDeriv, center = TRUE, scale = FALSE)
		sv <- svd(X)
		loadings <- sv$v                  # equals prcomp(...)$rotation
	}
	theOrder <- order(loadings[, 1])
	theCl = floor(seq(1, Nclusters+0.999, len= Nstrata))
	km = list(cluster = theCl)

	n_workers <- ceiling(config$num_threads / 2)

	strataMatrix = Matrix::sparseMatrix(
		i = seq_len(Nstrata),
		j = km$cluster,
		dims = c(Nstrata, max(km$cluster))
	)
	firstDerivBlock = firstDeriv %*% strataMatrix
	sparsityByBlockFirst = function(x, Sgamma0) {
		 whichAll = which(x != 0)-1L
  		whichRandom = as.vector(na.omit(match(whichAll, Sgamma0))-1L)
  		list(full = whichAll, random=whichRandom)
	}
	sparsityFirst = apply(firstDerivBlock, 2,sparsityByBlockFirst, Sgamma0)


	strataMatrixList = list(i=strataMatrix@i, j = as(strataMatrix, 'TsparseMatrix')@j, p=strataMatrix@p)

	strataListForHessian = mapply(
		function(x, y) {
			list(list(p = y$p[c(x, x+1)], i=y$i))
		},
		x = seq(1,len=ncol(strataMatrix)), 
		MoreArgs = list(y=strataMatrixList))

#	hessianDenseAll = hessianDenseLogical(parameters=x, data=dataNoMap,config=config, strata = res$groups$groups)
	# now find hessian for each block
	if(verbose) {
		cat("getting hessian by group...")
	}
	hessianByBlock = mapply(
		function(strata, config, ...) {
			hessianDenseLogical(..., config = modifyList(config, list(groups = strata)))
		},
		strata = strataListForHessian,
		MoreArgs = list(
			parameters=x, data=dataNoMap, 
			config = modifyList(config, list(verbose=FALSE, num_threads=1))),
		SIMPLIFY=FALSE#, mc.cores=config$num_threads
	)	
	hessianByBlock2 = lapply(hessianByBlock, function(xx) {
		xout = Matrix::Matrix(xx + t(xx), sparse=TRUE)
		xout@x = rep(1, length(xout@x))
		xout
	}
	)
	if(!all(unlist(lapply(hessianByBlock2, class)) == 'dsCMatrix')) {
		warning("some hessians not symmetric")
	} 
	if(!all(unlist(lapply(hessianByBlock2, function(xx) xx@uplo)) == 'U')){
		warning("some hessians not upper")
	} 

	#hessian for random part,
	if(verbose) {
		cat("getting Q hessian")
	}

	hessianQ = list(dense=hessianQdense(parameters=x, data=data, config=config))


# get non-zeros of tensor

	hessianQ$hessian = Matrix::forceSymmetric(Matrix::Matrix(hessianQ$dense + t(hessianQ$dense), 
		sparse=TRUE), uplo='U')
	hessianQT = as(as(hessianQ$hessian, 'generalMatrix'),'TsparseMatrix')
	hessianQns = as(hessianQ$hessian, 'generalMatrix')
	if(verbose) 
		cat("getting third sparsity...")

	ijkQ = hpolcc:::getThirdFromHessian(hessianQ$hessian)

	if(verbose) {
		cat("done\n")
	}
	if(any(ijkQ[,'Nunique'] == 3)) {
		warning("need to implement non-diagonal Q third")
	}
	fullHessian=do.call(rbind, lapply(hessianByBlock2, function(xx) cbind(xx@i,as(xx, "TsparseMatrix")@j)))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = hpolcc:::rowSortP(fullHessian) #t(apply(fullHessian, 1, sort))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = fullHessian[order(fullHessian[,2], fullHessian[,1]),]

	hessianNoQ = Matrix::sparseMatrix(i=apply(fullHessian,1,min), j=apply(fullHessian,1,max), 
		x=rep(1, nrow(fullHessian)), symmetric=TRUE, index1=FALSE)

	fullHessian = rbind(fullHessian, cbind(hessianQT@i, hessianQT@j))
	fullHessian = hpolcc:::rowSortP(fullHessian)
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	fullHessian = fullHessian[order(fullHessian[,2], fullHessian[,1]),]

	fullHessianMatrix = Matrix::sparseMatrix(
		i=fullHessian[,1], j=fullHessian[,2], 
		symmetric=TRUE, index1=FALSE,
		dims = c(Ntotal, Ntotal), repr='T')
	if(verbose) {
		cat("getting optimal pairs..")
	}

	fullList = hpolcc:::getOptimalPairs(
		fullHessianMatrix, 
		Sparams=Sparams, Sgamma1=Sgamma1,
		third=FALSE)

	if(verbose) {
		cat("done\n")
	}


	fullHessianPairs = paste(fullList$second$full$i, fullList$second$full$j, sep='_')
	fullHessianPairsNs = paste(fullList$second$nonSymmetric$i,fullList$second$nonSymmetric$j, sep='_')
	fullHessianPairsR = paste(fullList$second$random$i, fullList$second$random$j, sep='_')
	fullHessianPairsRNs = paste(fullList$second$randomNS$i, fullList$second$randomNS$j, sep='_')

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
	sparsityQ$nonSymmetric$match = try(match(
		paste(sparsityQ$nonSymmetric$i, sparsityQ$nonSymmetric$j, sep='_'), 
		fullHessianPairsNs
	)-1L)
	sparsityQ$random$match = try(match(
		paste(sparsityQ$random$i, sparsityQ$random$j, sep='_'), 
		fullHessianPairsR
	)-1L)
	sparsityQ$randomNS$match = try(match(
		paste(sparsityQ$randomNS$i, sparsityQ$randomNS$j, sep='_'), 
		fullHessianPairsRNs
	)-1L)


	# find full hessian sparsity
	# for each strata, get index in full hessian
	if(verbose) {
		cat("getting sparsity by block (this can take a while)...")				
	}
#				save(hessianByBlock2, Sparams, Sgamma1, fullHessianPairs, fullHessianPairsR, fullHessianPairsNs, fullHessianPairsRNs, file='todebug.Rdata')


	sparsityList = parallel::mcmapply(
		hpolcc:::getOptimalPairs,
		hessian = hessianByBlock2,
		grad = sparsityFirst,
		MoreArgs = list(Sparams = Sparams, Sgamma1=Sgamma1, 
			hessianPairs = fullHessianPairs,
			hessianPairsR = fullHessianPairsR, 
			hessianPairsNs = fullHessianPairsNs,
			hessianPairsRns = fullHessianPairsRNs
		), SIMPLIFY=FALSE, mc.cores=n_workers)

	if(verbose) {
		cat("done\n")
	}

	theTijk = lapply(sparsityList, function(xx) xx$third$ijk[,c('i','j','k')])
	theTijk = as.matrix(do.call(rbind, theTijk))
	theTijk = theTijk[!duplicated(theTijk), ]
	theTijk = rowSortP(theTijk)
	theTijk = theTijk[!duplicated(theTijk), ]
	colnames(theTijk) = c('i','j','k')
	theTijk = theTijk[
	order(theTijk[,'i'], theTijk[,'j'], theTijk[,'k']), 
	]
	fullList$third = list(ijk = theTijk)

	ijkChar = apply(theTijk, 1, paste, collapse='_')
	for(D in 1:length(sparsityList)) {
		ijkHere = rowSortP(as.matrix(sparsityList[[D]]$third$ijk[,c('i','j','k')]))
		sparsityList[[D]]$third$ijk$match =
		match(
			apply(ijkHere, 1, paste, collapse='_'),
			ijkChar
		) -1L
		if(any(is.na(sparsityList[[D]]$third$ijk$match))) {
			warning(D, "unmatched ijk")
		}
	}


	if('try-error' %in% class(sparsityList)) {
		sparsityList = list(
			hessianByBlock2=hessianByBlock2, Sparams= Sparams, Sgamma1=Sgamma1,
			hessianPairs = fullHessianPairs,
			hessianPairsR = fullHessianPairsR, 
			hessianPairsNs = fullHessianPairsNs,
			hessianPairsRns = fullHessianPairsRNs)
	}

	indexForThird = mapply(
		indexThirdTensor,
		i = seq(0L, len=Ntotal),
		MoreArgs = list(ijk = fullList$third$ijk),
		SIMPLIFY=FALSE
	)

	indexForDiag = mapply(
		indexThirdTensorDiag,
		i = seq(0L, len=Ntotal),
		MoreArgs = list(nonSymmetric = fullList$second$nonSymmetric),
		SIMPLIFY=FALSE
	)

	fullList$third$index = list(
		diag = indexForDiag,
		offDiag = indexForThird
	)

	result = list(
		group_sparsity = sparsityList,
		groups = strataMatrixList,
		sparsity = c(fullList, list(Q=sparsityQ)),
		firstDeriv = firstDeriv
	)

}