#' @export
adFun_groups = function(Ngroups, ATp) {
if(requireNamespace("RSpectra", quietly=TRUE) ) {
	loadings <- RSpectra::svds(ATp, k = 1)$v
} else {
	sv <- svd(ATp)
	loadings <- sv$v                  # equals prcomp(...)$rotation
}
theOrder <- order(loadings[, 1])
groupCut = seq(0, length(loadings)+1, len=Ngroups+1)

groupMat = Matrix::sparseMatrix(
	i = seq(0,len=length(loadings)),
	j = as.integer(cut(theOrder, groupCut))-1L,
	x=rep(1, length(loadings)), index1=FALSE)
groupMat
}
