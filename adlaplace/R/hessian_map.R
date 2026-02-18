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
	Sgamma0 = seq.int(Nbeta, length.out=Ngamma)
	Sgamma1 = Sgamma0+1L

	list_outer = list_inner = vector("list", Ngroups)
	for(D in seq_len(Ngroups)) {
		Nouter = length(sparsity_list[[D]]$col_hess)
		Ninner = length(sparsity_list[[D]]$col_hess_inner)

		list_outer[[D]] = data.table::data.table(
			shard = rep.int(D - 1L, Nouter),
			row = as.integer(sparsity_list[[D]]$row_hess),
			col = as.integer(sparsity_list[[D]]$col_hess),
			local = seq.int(0L, length.out = Nouter)
		)

		if (Ninner > 0L) {
			list_inner[[D]]  = data.table::data.table(
				shard = rep.int(D - 1L, Ninner),
				row = as.integer(sparsity_list[[D]]$row_hess_inner),
				col = as.integer(sparsity_list[[D]]$col_hess_inner),
				local = seq.int(0L, length.out = Ninner)
			)
		}
	}

	sparsityOuter = data.table::rbindlist(list_outer, use.names = TRUE)
	sparsityInner = data.table::rbindlist(list_inner, use.names = TRUE)

	sparsityOuter$cell = sparsityOuter$row + sparsityOuter$col * Nparams
	sparsityOuter$cellSparse = as.integer(factor(sparsityOuter$cell)) -1L


	sparsityInner$rowInner = match(sparsityInner$row, Sgamma0) -1L
	sparsityInner$colInner = match(sparsityInner$col, Sgamma0) -1L
	sparsityInner$cell = sparsityInner$rowInner + Ngamma*sparsityInner$colInner
	sparsityInner$cellSparse = as.integer(factor(sparsityInner$cell)) -1L


	outerMat = sparsityOuter[!duplicated(sparsityOuter$cell),]
	innerMat = sparsityInner[!duplicated(sparsityInner$cell),]

	hessian = Matrix::sparseMatrix(
		i=outerMat$row, j=outerMat$col,
		x = outerMat$cellSparse,
		symmetric=TRUE, index1=0, dims = rep(Nparams,2))
	# table(diff(hessian@x)) should be all ones

	hessian_inner = Matrix::sparseMatrix(
		i=innerMat$rowInner, j=innerMat$colInner,
		x = innerMat$cellSparse,
		symmetric=TRUE, index1=0, dims = rep(Ngamma,2))
	# table(diff(hessian_inner@x)) #should be all ones

	hessian_outer_map = Matrix::sparseMatrix(
		i = sparsityOuter$cellSparse, 
		j = sparsityOuter$shard,
		x = sparsityOuter$local,
		index1 = FALSE, 
		dims = c(length(unique(sparsityOuter$cellSparse)), Ngroups)
	)
	hessian_inner_map = Matrix::sparseMatrix(
		i = sparsityInner$cellSparse, 
		j = sparsityInner$shard,
		x = sparsityInner$local,
		index1 = FALSE, 
		dims = c(length(unique(sparsityInner$cellSparse)), Ngroups)
	)


	result_map = list(
		outer = list(
			p = as.integer(hessian_outer_map@p),
			local = as.integer(hessian_outer_map@x),
			global = as.integer(hessian_outer_map@i),
			dims = dim(hessian_outer_map)
		),

		inner = list(
			p = as.integer(hessian_inner_map@p),
			local = as.integer(hessian_inner_map@x),
			global = match(hessian_inner_map@i, as.integer(hessian_inner@x)) -1L,
			dims = c(length(hessian_inner@x), ncol(hessian_inner_map))
		)
	)

	list(
		hessian = list(outer = hessian, inner = hessian_inner),
		map = result_map
	)
}
