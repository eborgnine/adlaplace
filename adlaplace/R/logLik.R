

#' @export
logLik = function(x, data, config, 
	start_gamma = config$start_gamma, 	
	control = list(report.level=4, report.freq=1), 
	adPack, deriv=TRUE, 
	package = c(config$package, 'adlaplace')[1]
) {

	Nbeta = nrow(data$XTp)
	config_inner = config
	config_inner$beta = x[seq(1, len=Nbeta)]
	config_inner$theta = x[seq(Nbeta+1, len=length(x) - Nbeta)]
	Sgamma1 = seq(Nbeta+1, len=length(start_gamma))

	if(length(start_gamma) != nrow(data$ATp)) {
		warning("start_gamma is the wrong size")
	}

	if(config$verbose) {
		cat("logLik using package ", package, "for objective funcion\n")
	}

	if(missing(adPack)) {
		if(deriv) {
			adPack = getExportedValue(package, "getAdFun")(data, config, inner=FALSE)
		} else {
			adPack = NULL
		}
	}


	inner_res = getExportedValue(package, "inner_opt")(
		start_gamma,
		data=data, 
		config=config_inner, 
		control=control,
		adPackFull = adPack)

	names(inner_res$solution) = rownames(data$ATp)
	inner_res$parameters = x
	inner_res$full_parameters = c(config_inner$beta, inner_res$solution, config_inner$theta)
	try(names(inner_res$full_parameters) <- c(
		rownames(data$XTp), 
		rownames(data$ATp), 
		colnames(data$map)
	))


	result = c(
		list(
			inner = inner_res[grep("[pP]arameters|_full$", names(inner_res), invert=TRUE)],
			outer = inner_res[grep("_full$", names(inner_res))]
		),
		inner_res[grep("[pP]arameters|minusLogLik", names(inner_res))]
	)
	names(result$outer) = gsub("_outer$", "", names(result$outer))

	if(!deriv) {
		return(result)
	}

	result$deriv = logLikDeriv(
		x= result$full_parameters,
		inner_res = inner_res,
		config = config, 
		adPack = adPack
	)

	result$dLogLik = result$deriv$deriv$dL

return(result)
}

