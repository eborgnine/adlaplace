#' Find groups of observations with similar sparsity
#'
#' The function partitions columns into
#' \code{Ngroups} groups using quantiles of the first right singular vector,
#' optionally using \pkg{RSpectra} for efficiency when available.
#'
#' The resulting grouping is returned as a sparse matrix whose columns
#' correspond to groups and whose entries are the singular-vector loadings.
#' Groups are ordered from most heterogeneous to most homogeneous
#'
#' @param ATp A numeric matrix (or matrix-like object) whose columns are to
#'   be grouped. Typically this is a cross-product or design-related matrix.
#' @param Ngroups Integer giving the maximum number of groups to construct.
#'   The actual number of groups may be smaller if fewer distinct loadings
#'   are present.
#'
#' @details
#' If the \pkg{RSpectra} package is available, the leading singular vector
#' is computed using \code{RSpectra::svds}; otherwise, a full singular value
#' decomposition via \code{\link[base]{svd}} is used.
#'
#' Group boundaries are defined by empirical quantiles of the loadings.
#' Groups are subsequently reordered so that groups with larger within-group
#' variability appear first.
#'
#' @return
#' A sparse matrix of class \code{"dgCMatrix"} (from \pkg{Matrix}), with
#' one column per group and one row per column of \code{ATp}. Nonzero
#' entries correspond to singular-vector loadings.
#'
#'
#' @examples
#' set.seed(1)
#' A <- matrix(rnorm(100), 20, 5)
#' G <- adFun_groups(A, 3)
#' G
#'
#' @export
adFun_groups = function(ATp, elgm_matrix, Ngroups = ncol(ATp)) {

  if(!missing(elgm_matrix)) {
    ATp_t = methods::as(ATp, "TsparseMatrix")
    ATp_t = data.frame(row = ATp_t@j, gamma = ATp_t@i)
    elgm_matrix_t = methods::as(elgm_matrix, "TsparseMatrix")
    elgm_matrix_t = data.frame(row = elgm_matrix_t@i, strata=elgm_matrix_t@j)
    elgm_split = split(elgm_matrix_t, elgm_matrix_t$strata)
elgm_groups <- lapply(elgm_split, function(df, ATp_t) {
  Ahere <- ATp_t[ATp_t$row %in% df$row, ]
  data.frame(
    gamma = sort(unique(Ahere$gamma)),
    strata = df$strata[1]
  )
}, ATp_t = ATp_t)
  elgm_groups = do.call(rbind, elgm_groups)
  ATp = Matrix::sparseMatrix(i = elgm_groups$gamma, j=elgm_groups$strata,
    index1=FALSE, dims = c(nrow(ATp), ncol(elgm_matrix)))
  }

	if (inherits(ATp, "ngCMatrix")) {
 		 ATp <- methods::as(ATp, "dMatrix")  # numeric sparse, no dense copy
	}
	if(requireNamespace("RSpectra", quietly=TRUE) & (
    nrow(ATp) > 3
  ) ) {
		loadings <- RSpectra::svds(ATp, k = 1)$v[,1]
	} else {
		if(max(dim(ATp) > 1e5)) {
			warning("ATp matrix is very large, consider installing the RSpectra package")
		}
		sv <- svd(ATp)
		loadings <- sv$v[,1]                 # equals prcomp(...)$rotation
	}

	uniqueLoadings = unique(loadings)
	Ngroups = pmin(length(uniqueLoadings), Ngroups)
	groupCut = unique(stats::quantile(uniqueLoadings, seq(0,1,len=Ngroups+1)))
	groupCut[1] = groupCut[1]-1;groupCut[length(groupCut)] = groupCut[length(groupCut)]+1
	loadingsCut = cut(loadings, groupCut)
	loadingsCut2 = factor(loadingsCut, names(sort(table(loadingsCut), decreasing=TRUE)))
	theJ = as.integer(loadingsCut2)-1L

# order the groups from most heterogeneous to most homogeneous
	theSd = pmax(tapply(loadings, theJ, stats::sd),0, na.rm=TRUE)
	orderSd = order(theSd, decreasing=TRUE)-1L
	theJ2 = match(theJ, orderSd) -1L

	groupMat = Matrix::sparseMatrix(
		i = seq(0,len=length(loadings)),
		j = theJ2,
		x=	loadings, #rep(1L, length(loadings)), 
		index1=FALSE)
#plot(as(groupMat, 'TsparseMatrix')@j, groupMat@x)
	groupMat
}
