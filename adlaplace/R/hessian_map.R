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

	sparsity_list2 = vector("list", Ngroups)
	for(D in seq_len(Ngroups)) {
		if(config$verbose) {
			if(D %% 10000 == 0) {
				cat(" ", D, " ")
			}
		}
		Nouter = length(sparsity_list[[D]]$col_hess)
		Ninner = length(sparsity_list[[D]]$col_hess_inner)

		outer_dt = data.table::data.table(
			shard = rep.int(D - 1L, Nouter),
			inner = rep.int(0L, Nouter),
			row = as.integer(sparsity_list[[D]]$row_hess),
			col = as.integer(sparsity_list[[D]]$col_hess),
			local = seq.int(0L, length.out = Nouter)
		)

		if (Ninner > 0L) {
			inner_dt = data.table::data.table(
				shard = rep.int(D - 1L, Ninner),
				inner = rep.int(1L, Ninner),
				row = as.integer(sparsity_list[[D]]$row_hess_inner),
				col = as.integer(sparsity_list[[D]]$col_hess_inner),
				local = seq.int(0L, length.out = Ninner)
			)
			sparsity_list2[[D]] = data.table::rbindlist(list(outer_dt, inner_dt), use.names = TRUE)
		} else {
			sparsity_list2[[D]] = outer_dt
		}
	}
	if(config$verbose) {
		cat("done")
	}
	sparsityDf = data.table::rbindlist(sparsity_list2, use.names = TRUE)

	fullMat = unique(sparsityDf[, .(inner, row, col)])

	outerMat = fullMat[inner == 0L]
	hessian = Matrix::sparseMatrix(
		i=outerMat$row, j=outerMat$col,
		x = rep(-1, nrow(outerMat)),
		symmetric=TRUE, index1=0, dims = rep(Nparams,2))
	hessian@x = seq(0, length.out = length(hessian@x))

	hessianT = as(hessian, "TsparseMatrix")
	indexMat_outer = data.table::data.table(
		index = as.integer(hessianT@x),
		row = as.integer(hessianT@i),
		col = as.integer(hessianT@j),
		inner = 0L
	)


	hessian_inner = hessian[Sgamma1, Sgamma1]
	indexMat_inner = indexMat_outer[row %in% Sgamma0 & col %in% Sgamma0]
	indexMat_inner[, inner := 1L]

	indexMat = data.table::rbindlist(list(indexMat_outer, indexMat_inner), use.names = TRUE)
	if(config$verbose) {
		cat(" merge ")
	}

	data.table::setkey(indexMat, inner, row, col)
	data.table::setkey(sparsityDf, inner, row, col)
	indexDf <- data.table::merge.data.table(
		indexMat, sparsityDf,
		all = TRUE,
		sort = FALSE
	)
	if(config$verbose) {
		cat("done.")
	}
	if (indexDf[, any(is.na(index) | is.na(shard))]) {
		warning("problem merging hessian indices")
	}


	build_shard_map = function(xx, Nglobal, Ngroups) {
		xx = xx[!is.na(index) & !is.na(shard) & !is.na(local)]
		data.table::setorderv(xx, c("shard", "local"))
		theTable = tabulate(xx$shard + 1L, nbins = Ngroups)
		Matrix::sparseMatrix(
			p = c(0L, cumsum(as.integer(theTable))),
			i = as.integer(xx$index),
			x = as.integer(xx$local),
			index1 = FALSE,
			dims = c(Nglobal, Ngroups)
		)
	}

	outerDf = indexDf[inner == 0L]
	innerDf = indexDf[inner == 1L]
	hessian_outer_map = build_shard_map(outerDf, length(hessian@x), Ngroups)
	hessian_inner_map = build_shard_map(innerDf, length(hessian@x), Ngroups)

	if(config$verbose) {
		cat("-")
	}

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
			global = match(hessian_inner_map@i, hessian_inner@x) -1L,
			dims = c(length(hessian_inner@x), ncol(hessian_inner_map))
		)
	)


	result_map
	list(
		hessian = list(outer = hessian, inner = hessian_inner),
		map = result_map
	)
}
