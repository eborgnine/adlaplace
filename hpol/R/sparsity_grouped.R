#' @export

sparsity_grouped = function(x, data, config, verbose=FALSE) {

# library('hpolcc');x = res$parameters_and_gamma;data=res$tmb_data;config=res$config;verbose=TRUE

	Nstrata = ncol(data$cc_matrixTp)
	Sstrata = seq(0, len=Nstrata)

	Nclusters = config$Nclusters
	if(!length(Nclusters)) Nclusters  <- config$num_threads*2                                   

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

			n_workers <- ceiling(config$num_threads / 2)
			cl <- parallel::makeCluster(n_workers, type = "PSOCK")
#			on.exit(parallel::stopCluster(cl), add = TRUE)

		if(Nclusters == 1) {
			km = list(cluster = rep(1, Nstrata))		
		} else {

			tFirst <- t(firstDeriv)                 # do this once in the master
			centers <- ceiling(1.1 * Nclusters)
			nstart  <- max(5L, ceiling(2 * Nclusters / config$num_threads))

			parallel::clusterExport(cl, varlist = c("tFirst", "centers", "nstart"), envir = environment())
			parallel::clusterEvalQ(cl, { gc(); NULL })       # optional hygiene
			parallel::clusterSetRNGStream(cl, 123)           # reproducible, different stream per worker
			seeds <- seq_len(n_workers)
			kmMC = parallel::parLapply(cl, seeds, function(seed) {
				set.seed(seed)
				stats::kmeans(tFirst, centers = centers, iter.max = 1000,
					nstart = nstart, algorithm = "Hartigan-Wong")
			})
			parallel::stopCluster(cl)
			
			km <- kmMC[[ which.min(sapply(kmMC, `[[`, "tot.withinss")) ]]

#	km <- try(kmeans(t(firstDeriv), centers = ceiling(1.1*Nclusters), 
#		iter.max=25000, nstart=10*Nclusters, algorithm='Hartigan-Wong'))

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
						maxClustersInd = order(theTable, decreasing=TRUE)[1:2]
						maxClustersSize = theTable[maxClustersInd]
						if( maxClustersSize[1] > 1.5*maxClustersSize[2]) {
							whichInBiggest = which(km$cluster == maxClustersInd[1])
							firstDeriv2 = firstDeriv[,whichInBiggest]
							pr = prcomp(firstDeriv2)
							theOrder = order(pr$rotation[,1])
							newk = cut(theOrder, seq(0, length(theOrder)+1, 
								len=ceiling(maxClustersSize[1]/maxClustersSize[2])))
							km$cluster[whichInBiggest] = -as.numeric(newk)
							km$cluster = as.integer(factor(km$cluster))
						}
					}
				} else {
					pr = prcomp(firstDeriv)
					theOrder = order(pr$rotation[,1])
					theCl = floor(seq(1, Nclusters+0.999, len= Nstrata))
					km = list(cluster = theCl)
				}
			} # not one cluster

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
				hessianByBlock = mapply(function(strata, config, ...) {
					hessianDenseLogical(..., config = c(config, list(groups = strata)))
				},
				strata = strataListForHessian,
				MoreArgs = list(
					parameters=x, data=dataNoMap, 
					config = c(config[setdiff(names(config), 'num_threads')], 
						list(num_threads=1))),
				SIMPLIFY=FALSE#, mc.cores=config$num_threads
			)
			hessianByBlock2 = lapply(hessianByBlock, Matrix::Matrix, sparse=TRUE)

	#hessian for random part,
			if(verbose) cat("getting Q hessian")
				hessianQ = list(dense=hessianQdense(parameters=x, data=data, config=config))


# get non-zeros of tensor

			hessianQ$hessian = Matrix::forceSymmetric(Matrix::Matrix(hessianQ$dense, sparse=TRUE), uplo='U')
			hessianQT = as(as(hessianQ$hessian, 'generalMatrix'),'TsparseMatrix')
			hessianQns = as(hessianQ$hessian, 'generalMatrix')
			if(verbose) cat("getting third sparsity...")
				ijkQ = hpolcc:::getThirdFromHessian(hessianQ$hessian)
			if(verbose) cat("done\n")
				if(any(ijkQ[,'Nunique'] == 3)) warning("need to implement non-diagonal Q third")

					fullHessian=do.call(rbind, lapply(hessianByBlock2, function(xx) cbind(xx@i,as(xx, "TsparseMatrix")@j)))
				fullHessian = fullHessian[!duplicated(fullHessian), ]
				fullHessian = rowSortP(fullHessian) #t(apply(fullHessian, 1, sort))
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

				fullList = hpolcc:::getOptimalPairs(fullHessianMatrix, Sparams=Sparams, Sgamma1=Sgamma1)

				if(verbose) cat("done\n")


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
if(verbose) cat("getting sparsity by block...")				
				save(hessianByBlock2, Sparams, Sgamma1, fullHessianPairs, fullHessianPairsR, fullHessianPairsNs, fullHessianPairsRNs, file='todebug.Rdata')

if(FALSE) {
				parallel::clusterExport(cl, c( "hessianByBlock2", "Sparams", "Sgamma1",
                    "fullHessianPairs", "fullHessianPairsR",
                    "fullHessianPairsNs", "fullHessianPairsRNs"), envir=environment())
parallel::clusterEvalQ(cl, { library(hpolcc); NULL })



					sparsityList = parallel::parLapply(cl, seq_along(hessianByBlock2), function(i) {
					  hpolcc:::getOptimalPairs(
						hessian = hessianByBlock2[[i]],
						Sparams = Sparams, Sgamma1=Sgamma1, 
							hessianPairs = fullHessianPairs,
							hessianPairsR = fullHessianPairsR, 
							hessianPairsNs = fullHessianPairsNs,
							hessianPairsRns = fullHessianPairsRNs
						)
					}
					)
}
				
				sparsityList = mapply(
				  hpolcc:::getOptimalPairs,
				    hessian = hessianByBlock2,
				    MoreArgs = list(Sparams = Sparams, Sgamma1=Sgamma1, 
				    hessianPairs = fullHessianPairs,
				    hessianPairsR = fullHessianPairsR, 
				    hessianPairsNs = fullHessianPairsNs,
				    hessianPairsRns = fullHessianPairsRNs
				  ), SIMPLIFY=FALSE)
				
				if(verbose) cat("done\n")



					if('try-error' %in% class(sparsityList)) sparsityList = list(hessianByBlock2=hessianByBlock2, Sparams= Sparams, Sgamma1=Sgamma1,
							hessianPairs = fullHessianPairs,
							hessianPairsR = fullHessianPairsR, 
							hessianPairsNs = fullHessianPairsNs,
							hessianPairsRns = fullHessianPairsRNs)

					result = list(
						group_sparsity = sparsityList,
						groups = strataMatrixList,
						sparsity = c(fullList, list(Q=sparsityQ)),
						firstDeriv = firstDeriv
					)

			}