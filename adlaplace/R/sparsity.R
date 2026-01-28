#' Build grouped Hessian sparsity structures
#'
#' Constructs Hessian sparsity “templates” (global/outer and inner-\code{gamma})
#' and then derives per-group sparsity objects aligned to those templates.
#'
#'
#' @param data A list containing model matrices and metadata required by
#'   \code{sparsity()} and by downstream code. (Backend-dependent.)
#' @param config A list of configuration options. This function expects at least
#'   \code{beta}, \code{start_gamma}, \code{theta}, and optionally \code{num_threads}.
#' @param sparsity_list Optional list describing sparsity per group. Each element
#'   is expected to contain integer vectors \code{row} and \code{col} giving
#'   Hessian nonzero coordinates (0-based indices) for that group.
#'
#'
#' @return A list with components:
#' \describe{
#'   \item{hessian_outer}{Symmetric sparse matrix template (\code{dgCMatrix}) for the full parameter Hessian (upper triangle).}
#'   \item{hessianL_outer}{Lower-triangular variant of \code{hessian_outer}.}
#'   \item{group_outer}{List of per-group sparsity objects aligned to \code{hessian_outer}.}
#'   \item{hessian_inner}{Sparse matrix template for the \code{gamma} block of the Hessian.}
#'   \item{hessianL_inner}{Lower-triangular variant of \code{hessian_inner}.}
#'   \item{group_inner}{List of per-group sparsity objects aligned to \code{hessian_inner}.}
#' }
#'
#' @examples
#' \dontrun{
#' # sp <- group_sparsity(data, config)
#' # str(sp$hessian_outer)
#' # length(sp$group_outer)
#' }
#'


#' @export
group_sparsity = function(data, config, sparsity_list) {

	if(missing(sparsity_list)) {
		package = c(config$package, 'adlaplace')[1]
		sparsity_list = 
			getExportedValue(package, "sparsity")(data, config)
	}

	Sgamma0 = seq.int(length(config$beta), length.out=length(config$start_gamma))
	Sgamma1 = Sgamma0+1L
	Nparams = length(config$beta) + length(config$start_gamma) + length(config$theta)


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
	hessianTemplateT = methods::as(hessianTemplate, 'TsparseMatrix')

	hessianTemplateInner = hessianTemplate[Sgamma1, Sgamma1]
	hessianTemplateInner@x = seq(0L, len=length(hessianTemplateInner@x))
	hessianTemplateInnerL = Matrix::t(hessianTemplateInner)
	hessianTemplateInnerT = methods::as(hessianTemplateInner, 'TsparseMatrix')


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


sparsity_by_group = function(xx, template, dims, Sgamma0) {
	allIJ = cbind(i=xx$row, j=xx$col)
	if(!missing(Sgamma0)) {
		allIJ = cbind(
			i=match(allIJ[,'i'], Sgamma0)-1L,
			j=match(allIJ[,'j'], Sgamma0)-1L
		)
		allIJ = allIJ[which(!(is.na(allIJ[,'i']) | is.na(allIJ[,'j']))), ,drop=FALSE]
		xx$grad = as.vector(stats::na.omit(match(xx$grad, Sgamma0) -1L))
	}
	allIJ = allIJ[!duplicated(allIJ), ,drop=FALSE]
	allIJ = allIJ[which(allIJ[,'j'] >= allIJ[,'i']),,drop=FALSE]
	allIJ = allIJ[order(allIJ[,'j'], allIJ[,'i']),,drop=FALSE]

	hessianTemplateHere = Matrix::sparseMatrix(
		i=allIJ[,'i'], j=allIJ[,'j'],
		symmetric=TRUE, index1=FALSE,
		dims = dims
	)
	hessianTemplateHereT = methods::as(hessianTemplateHere, 'TsparseMatrix')
	ijHere = paste0(hessianTemplateHereT@i, rep('_', length(hessianTemplateHereT@i)),
		hessianTemplateHereT@j)

	matchHere = as.vector(stats::na.omit(match(ijHere, template)))-1L


	list(
			grad = xx$grad,
			i=hessianTemplateHereT@i, 
			j=hessianTemplateHereT@j,
			match = matchHere,
		p = hessianTemplateHere@p
	)
}
