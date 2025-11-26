
sparsity_by_group = function(xx, template, dims) {
	allIJ = cbind(i=xx$row, j=xx$col)
	allIJ = allIJ[!duplicated(allIJ), ]
	allIJ = allIJ[allIJ[,'j'] >= allIJ[,'i'],]
	allIJ = allIJ[order(allIJ[,'j'], allIJ[,'i']),]

	hessianTemplateHere = Matrix::sparseMatrix(
		i=allIJ[,'i'], j=allIJ[,'j'],
		symmetric=TRUE, index1=FALSE,
		dims = dims
	)
	hessianTemplateHereT = as(hessianTemplateHere, 'TsparseMatrix')
	ijHere = paste0(hessianTemplateHereT@i, '_', hessianTemplateHereT@j)

	matchHere = match(ijHere, template)-1L
	list(
			grad = xx$grad,
			i=hessianTemplateHereT@i, 
			j=hessianTemplateHereT@j,
			match = matchHere,
		p = hessianTemplateHere@p
	)
}

#' @export
group_sparsity = function(adFun, config) {

	sparsity_list = sparsity(adFun, config$start_gamma)
	allRow = lapply(sparsity_list, '[[', 'row')
	allCol = lapply(sparsity_list, '[[', 'col')
	allIJ = cbind(i=unlist(allRow), j=unlist(allCol))
	allIJ = allIJ[!duplicated(allIJ), ]
	allIJ = allIJ[allIJ[,'j'] >= allIJ[,'i'],]
	allIJ = allIJ[order(allIJ[,'j'], allIJ[,'i']),]
	hessianTemplate = Matrix::sparseMatrix(
		i=allIJ[,'i'], j=allIJ[,'j'],
		symmetric=TRUE, index1=FALSE,
		dims = rep(length(config$start_gamma),2)
	)
	hessianTemplateT = as(hessianTemplate, 'TsparseMatrix')
	templateIJ = paste0(hessianTemplateT@i, '_', hessianTemplateT@j)
	sparsity_list2 = parallel::mclapply(sparsity_list, 
		sparsity_by_group, 
		template = templateIJ,
		dims = dim(hessianTemplate),
		mc.cores = 
			max(c(config$num_threads, 1), na.rm=TRUE)
	)
	return(list(hessian = hessianTemplate, group = sparsity_list2))
}