# map estimator
#' @export
map_estimator = function(
	parameters, 
	data, config,
	control=list()
) {




	wrappers <- make_trustoptim_wrappers(
		data = data,
		config = config[setdiff(names(config), 'sparsity')]
	)

	if(FALSE) {
		stuff = objectiveFunctionC(
			parameters = parameters, data=data, 
			config=config[setdiff(names(config), 'sparsity')])
		wrappers$fn(parameters)
		wrappers$fn(parameters+0.001)

		stuff1 = wrappers$gr(parameters)
		stuff2 = wrappers$hs(parameters)
		stuff3 = wrappers$hs(parameters)
	}

	if(!identical(config$transform_theta, TRUE)) {
		warning("transform_theta FALSE not implemented yet")
	}

  # inner opt
	result <- trustOptim::trust.optim(
		x = parameters,
		fn = wrappers$fn,
		gr = wrappers$gr,
		hs = wrappers$hs,
		method = "Sparse",
		control = control
	)
	result$map = list(
		parameters = list(
			beta = result$solution[1:nrow(data$XTp)],
			theta = result$solution[-seq(1, nrow(data$XTp) + nrow(data$ATp))]
		),
		random = result$solution[seq(
			from=nrow(data$XTp)+1, len=nrow(data$ATp)
		)]
	)

	result
}