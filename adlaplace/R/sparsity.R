
sparsity_by_group = function(xx, template, dims, Sgamma0) {
	allIJ = cbind(i=xx$row, j=xx$col)
	if(!missing(Sgamma0)) {
		allIJ = cbind(
			i=match(allIJ[,'i'], Sgamma0)-1L,
			j=match(allIJ[,'j'], Sgamma0)-1L
		)
		allIJ = allIJ[which(!(is.na(allIJ[,'i']) | is.na(allIJ[,'j']))), ]
		xx$grad = as.vector(na.omit(match(xx$grad, Sgamma0) -1L))
	}
	allIJ = allIJ[!duplicated(allIJ), ,drop=FALSE]
	allIJ = allIJ[which(allIJ[,'j'] >= allIJ[,'i']),,drop=FALSE]
	allIJ = allIJ[order(allIJ[,'j'], allIJ[,'i']),,drop=FALSE]

	hessianTemplateHere = Matrix::sparseMatrix(
		i=allIJ[,'i'], j=allIJ[,'j'],
		symmetric=TRUE, index1=FALSE,
		dims = dims
	)
	hessianTemplateHereT = as(hessianTemplateHere, 'TsparseMatrix')
	ijHere = paste0(hessianTemplateHereT@i, rep('_', length(hessianTemplateHereT@i)),
		hessianTemplateHereT@j)

	matchHere = as.vector(na.omit(match(ijHere, template)))-1L


	list(
			grad = xx$grad,
			i=hessianTemplateHereT@i, 
			j=hessianTemplateHereT@j,
			match = matchHere,
		p = hessianTemplateHere@p
	)
}

#' @export
group_sparsity = function(data, config, adFun) {

	Sgamma0 = seq.int(length(config$beta), length.out=length(config$start_gamma))
	Sgamma1 = Sgamma0+1L
	Nparams = length(config$beta) + length(config$start_gamma) + length(config$theta)

	if(missing(adFun)) {
		adFun = adlaplace::getAdFun(data, config, inner=FALSE)
	}

	sparsity_list = adlaplace:::sparsity(
		adFun, 
		c(config$beta, config$start_gamma, config$theta),
		verbose=pmax(FALSE,config$verbose, na.rm=TRUE))
	allRow = lapply(sparsity_list, '[[', 'row')
	allCol = lapply(sparsity_list, '[[', 'col')
	allIJ = cbind(i=unlist(allRow), j=unlist(allCol))
	allIJ = allIJ[!duplicated(allIJ), ]
	allIJ = allIJ[allIJ[,'j'] >= allIJ[,'i'],]
	allIJ = allIJ[order(allIJ[,'j'], allIJ[,'i']),]
	hessianTemplate = Matrix::sparseMatrix(
		i=allIJ[,'i'], j=allIJ[,'j'],
		x = seq(0L, len=nrow(allIJ)),
		symmetric=TRUE, index1=FALSE,
		dims = rep(Nparams,2)
	)
	hessianTemplateL = Matrix::t(hessianTemplate) # lower triangle version
	hessianTemplateT = as(hessianTemplate, 'TsparseMatrix')

	hessianTemplateInner = hessianTemplate[Sgamma1, Sgamma1]
	hessianTemplateInner@x = seq(0L, len=length(hessianTemplateInner@x))
	hessianTemplateInnerL = Matrix::t(hessianTemplateInner)
	hessianTemplateInnerT = as(hessianTemplateInner, 'TsparseMatrix')


	templateIJ = paste0(hessianTemplateT@i, '_', hessianTemplateT@j)
	templateIJinner = paste0(hessianTemplateInnerT@i, '_', hessianTemplateInnerT@j)

	sparsity_list_outer = parallel::mclapply(
		sparsity_list, 
		sparsity_by_group, 
		template = templateIJ,
		dims = dim(hessianTemplate),
		mc.cores = max(c(config$num_threads, 1), na.rm=TRUE)
	)
	sparsity_list_inner = parallel::mclapply(
		sparsity_list, 
		sparsity_by_group, 
		template = templateIJinner,
		dims = dim(hessianTemplateInner),
		Sgamma0 = Sgamma0,
		mc.cores = 
			max(c(config$num_threads, 1), na.rm=TRUE)
	)
	return(list(
		hessian_outer = hessianTemplate, hessianL_outer = hessianTemplateL, group_outer = sparsity_list_outer,
		hessian_inner = hessianTemplateInner, hessianL_inner = hessianTemplateInnerL, group_inner = sparsity_list_inner
	))
}