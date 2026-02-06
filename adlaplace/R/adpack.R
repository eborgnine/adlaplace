

getAdpackHandle = function(data, config) {
	theAdFun = getAdFun(data, config)
	theSparsity = extractSparsity(theAdFun)
	# this is the part that can't be done in cpp easily
	theHessianMap = hessianMap(config, theSparsity)
	getAdHandle(theAdFun, theHessianMap)


}