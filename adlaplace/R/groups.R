#' @export
adFun_groups = function(Ngroups, ATp) {
if(requireNamespace("RSpectra", quietly=TRUE) ) {
	loadings <- RSpectra::svds(ATp, k = 1)$v[,1]
} else {
	sv <- svd(ATp)
	loadings <- sv$v[,1]                 # equals prcomp(...)$rotation
}


#km = kmeans(loadings, centers=Ngroups)

uniqueLoadings = unique(loadings)
Ngroups = pmin(length(uniqueLoadings), Ngroups)
groupCut = unique(quantile(uniqueLoadings, seq(0,1,len=Ngroups+1)))
groupCut[1] = groupCut[1]-1;groupCut[length(groupCut)] = groupCut[length(groupCut)]+1
theJ = as.integer(cut(loadings, groupCut))-1L

# order the groups from most heterogeneous to most homogeneous
theSd = pmax(tapply(loadings, theJ, sd),0, na.rm=TRUE)
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
