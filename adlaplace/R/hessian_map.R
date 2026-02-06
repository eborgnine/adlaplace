#' Build grouped Hessian sparsity structures
#'
#' Constructs Hessian sparsity “templates” (global/outer and inner-\code{gamma})
#' and then derives per-group sparsity objects aligned to those templates.
#'
#'
#' @param config A list of configuration options. This function expects at least
#'   \code{beta}, \code{gamma}, \code{theta}.
#' @param sparsity_list List describing sparsity obtained from \link{getAdFun}.
#'
#'
#' @return A list with components \code{hessian, hessian_inner}, and \code{match} to be added to \code{config$sparsity}
#' \code{match$hessian} is a sparse matrix, entry A of \code{hessian@x} for shard \code{B} is \code{match$hessian[A,B]}
#'


#' @export
hessianMap = function(sparsity_list, config) {

	Nbeta = length(config$beta)
	Ngamma = length(config$gamma)
	Ntheta = length(config$theta)

	Ngroups = length(sparsity_list)
	Nparams = Nbeta + Ngamma + Ntheta
	Sgamma0 = seq.int(Nbeta length.out=Ngamma)
	Sgamma1 = Sgamma0+1L

	sparsity_list2 = list()
	for(D in 1:Ngroups) {
		Nhere = lapply(sparsity_list[[D]], length)
		repZero = lapply(Nhere, rep, x=0L)
		sparsity_list2[[D]] = cbind(
			shard = rep(D-1L, Nhere$col_hess + Nhere$col_hess_inner),
			rbind(
				cbind(
					inner=repZero$col_hess, 
					row=sparsity_list[[D]]$row_hess, 
					col=sparsity_list[[D]]$col_hess,
					local = seq.int(0L, length.out=Nhere$row_hess)),
				cbind(
					inner=1L + repZero$col_hess_inner, 
					row=sparsity_list[[D]]$row_hess_inner, 
					col=sparsity_list[[D]]$col_hess_inner,
					local = seq.int(0L, length.out=Nhere$row_hess_inner))
			)
		)
	}
	sparsityDf = do.call(rbind, sparsity_list2)

	fullMat = sparsityDf[, c('inner','row','col')] 
	fullMat = fullMat[!duplicated(fullMat),]

	outerMat = fullMat[fullMat[,'inner'] == 0, ]
	hessian = Matrix::sparseMatrix(
		i=outerMat[,'row'], j=outerMat[,'col'],
		x = rep(-1, nrow(outerMat)),
		symmetric=TRUE, index1=0, dims = rep(Nparams,2))
	hessian@x = seq(0, length.out = length(hessian@x))
	hessian_inner = hessian[Sgamma1, Sgamma1]

	hessianT = as(hessian, "TsparseMatrix")
	indexMat = data.frame(index=as.integer(hessianT@x), row=hessianT@i, col=hessianT@j)
	indexMat$index_inner = match(indexMat$index, as.integer(hessian_inner@x)) -1L

	indexDf = merge(indexMat, sparsityDf, all=TRUE)
	if(any(is.na(unlist(indexDf[c('index','shard')])))) {
		warning("problem merging hessian indices")
	}
	whichInner = which(indexDf$inner == 1)
	indexDf[whichInner, 'index'] = 	indexDf[whichInner, 'index_inner']

	indexDf$shardFac = factor(indexDf$shard, seq.int(0L, length.out=Ngroups))
	indexDf$innerFac = factor(indexDf$inner, levels = c(0L, 1L), labels = c('outer','inner'))
	indexDf = indexDf[order(indexDf$innerFac, indexDf$shard, indexDf$local), ]
#		c('index','innerFac','shardFac','local')

	indexHessianSplit = split(indexDf, indexDf$innerFac)
	hessianSparsity = lapply(
		indexHessianSplit, function(xx) {
			theTable = as.integer(table(xx$shardFac))
			Matrix::sparseMatrix(
				p = c(0L, cumsum(theTable)), 
				i = xx$index, 
				x = xx$local,
				index1=FALSE,
				dims = c(length(hessian@x), Ngroups))
		}
	)
#	names(hessianSparsity) = gsub("^outer$", "hessian_outer",names(hessianSparsity))
#	names(hessianSparsity) = gsub("^inner$", "hessian_inner",names(hessianSparsity))
	# D=2;bob = Matrix::sparseMatrix(i = hessianSparsity[[D]]$local, p=hessianSparsity[[D]]$p, x = as.numeric(hessianSparsity[[D]]$index), index1=0)

	result_map = lapply(hessian_sparsity, function(xx) {
		list(
			p = as.integer(xx@p),
			local = as.integer(xx@i),
			global = as.integer(xx@x),
			dims = dim(xx)
		)
	})

	list(
		hessian = list(outer = hessian, inner = hessian_inner),
		map = result_map
	)
}


