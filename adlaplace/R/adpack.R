
getBackendContext = function(
	adPack,
	map
) {


}

adFunList = getAdFun_api(data, config)
config$sparsity = adlaplace::group_sparsity(config, sparsity_list = adFunList$sparsity)

backendContext = getBackendContext(adFunList$adPack, config$sparsity$hessian_map)