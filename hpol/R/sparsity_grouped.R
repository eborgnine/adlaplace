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

config$verbose=TRUE
	firstDeriv = mapply(logLikNoQStrata,
		strata = singleStrataList,
		MoreArgs = list(
			parameters=x, data=data, sparsity = list(),
			config = c(config[
				setdiff(names(config), c('num_Rthreads', 'dense', 'maxDeriv'))
				], list(dense=TRUE, maxDeriv=1, num_threads=1))),
		SIMPLIFY=FALSE
	)
	firstDeriv = do.call(cbind, lapply(firstDeriv, function(xx) xx$grad))
	firstDeriv = firstDeriv > 0

		km <- try(kmeans(t(firstDeriv), centers = ceiling(1.5*Nclusters), iter.max=15000,nstart=5*Nclusters, algorithm='Hartigan-Wong'))
		if(!any(class('km') == 'try-error')) {
			mergeThreshold = floor(Nstrata/Nclusters)
			NtoMerge = max(which(cumsum(sort(km$size)) < mergeThreshold))
			mergeSize = sort(km$size)[NtoMerge]
			whichToMerge = which(km$size <= mergeSize)

			km$cluster2 = km$cluster
			km$cluster2[km$cluster %in% whichToMerge] = min(whichToMerge)
			km$cluster = as.integer(factor(km$cluster2))
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

	strataMatrixList = list(i=strataMatrix@i, p=strataMatrix@p)

	strataListForHessian = mapply(
		function(x, y) {
			list(list(p = y$p[c(x, x+1)], i=y$i))
		},
		x = seq(1,len=ncol(strataMatrix)), 
		MoreArgs = list(y=strataMatrixList))


	# now find hessian for each block
	hessianByBlock = parallel::mcmapply(logLikNoQStrata,
		strata = strataListForHessian,
		MoreArgs = list(
		parameters=x, data=data, sparsity = list(),
		config = c(config[setdiff(names(config), c('dense','num_threads', 'maxDeriv'))], list(maxDeriv=2, num_threads=1, dense=TRUE))),
		SIMPLIFY=FALSE, mc.cores=config$num_threads
	)
	hessianByBlock2 = lapply(hessianByBlock, function(xx) Matrix::Matrix(xx$hessian, sparse=TRUE))

	fullHessian=do.call(rbind, lapply(hessianByBlock2, function(xx) cbind(xx@i,as(xx, "TsparseMatrix")@j)))
	fullHessian = fullHessian[!duplicated(fullHessian), ]
	colnames(fullHessian) = c('i','j')
	fullHessian = fullHessian[fullHessian[,'i'] <= fullHessian[,'j'],]
	fullHessian = fullHessian[order(fullHessian[,'j'], fullHessian[,'i']),]


	# find full hessian sparsity
	# for each strata, get index in full hessian

	sparsityList = parallel::mcmapply(getOptimalPairs,
		hessian = hessianByBlock2,
		MoreArgs = list(Sparams = Sparams, Sgamma1=Sgamma1, fullHessian = fullHessian), 
		SIMPLIFY=FALSE, mc.cores=config$num_threads)


	# try to get the hessian
	bob = logLikNoQStrata(
		parameters=x, data=data, config = c(
			config[setdiff(names(config), c('sparsity','maxDeriv'))],
			list(maxDeriv=2)), 
		strata = strataMatrixList,
		sparsity = sparsityList
	)

}