#' Build grouped Hessian sparsity structures
#'
#' Constructs Hessian sparsity “templates” (global/outer and inner-\code{gamma})
#' and then derives per-group sparsity objects aligned to those templates.
#'
#'
#' @param sparsity_list List describing sparsity obtained from \link{getAdFun}.
#' @param Nbeta Number of fixed-effect parameters.
#' @param Ngamma Number of random-effect parameters.
#' @param Ntheta Number of theta parameters.
#'
#'
#' @return A list with components:
#' \describe{
#'   \item{hessian}{A list with \code{outer} and \code{inner} symmetric sparse
#'   template matrices (\code{dsCMatrix}). Their \code{x} slots encode
#'   contiguous global nonzero ids.}
#'   \item{map}{A list with \code{outer} and \code{inner} maps. Each map
#'   contains integer vectors \code{p}, \code{local}, \code{global}, and
#'   \code{dims} describing how each group-local Hessian nonzero maps to a
#'   global template nonzero.}
#' }
#'


#' @export
hessianMap = function(sparsity_list, Nbeta, Ngamma, Ntheta) {

	Ngroups = length(sparsity_list)
	Nparams = Nbeta + Ngamma + Ntheta
	Sgamma0 = seq.int(Nbeta, length.out=Ngamma)

	list_outer = list_inner = vector("list", Ngroups)
	make_sparse_df <- function(shard_id, row, col) {
		n <- length(col)
		data.frame(
			shard = rep.int(as.integer(shard_id), n),
			row = as.integer(row),
			col = as.integer(col),
			local = seq.int(0L, length.out = n),
			stringsAsFactors = FALSE
		)
	}
	bind_sparse_df <- function(lst) {
		lst <- lst[lengths(lst) > 0L]
		if (length(lst) == 0L) {
			return(data.frame(
				shard = integer(0),
				row = integer(0),
				col = integer(0),
				local = integer(0),
				stringsAsFactors = FALSE
			))
		}
		do.call(rbind, lst)
	}
	for(D in seq_len(Ngroups)) {
		list_outer[[D]] = make_sparse_df(
			shard_id = D - 1L,
			row = sparsity_list[[D]]$row_hess,
			col = sparsity_list[[D]]$col_hess
		)

		list_inner[[D]]  = make_sparse_df(
			shard_id = D - 1L,
			row = sparsity_list[[D]]$row_hess_inner,
			col = sparsity_list[[D]]$col_hess_inner
		)
	}

	sparsityOuter = bind_sparse_df(list_outer)
	sparsityInner = bind_sparse_df(list_inner)

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
		x = outerMat$cellSparse, repr='C',
		symmetric=TRUE, index1=0, dims = rep(Nparams,2))
	# table(diff(hessian@x)) should be all ones

	hessian_inner = Matrix::sparseMatrix(
		i=innerMat$rowInner, j=innerMat$colInner,
		x = innerMat$cellSparse, repr='C',
		symmetric=TRUE, index1=0, dims = rep(Ngamma,2))
	# table(diff(hessian_inner@x)) #should be all ones

	hessian_outer_map = Matrix::sparseMatrix(
		i = sparsityOuter$cellSparse, 
		j = sparsityOuter$shard,
		x = sparsityOuter$local,
		index1 = FALSE, repr='C',
		dims = c(length(unique(sparsityOuter$cellSparse)), Ngroups)
	)
	hessian_inner_map = Matrix::sparseMatrix(
		i = sparsityInner$cellSparse, 
		j = sparsityInner$shard,
		x = sparsityInner$local,
		index1 = FALSE, repr='C',
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
