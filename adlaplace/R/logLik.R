#' @export
logLik = function(x, data, config, 
	start_gamma = config$start_gamma, 	
	control = list(report.level=4, report.freq=1), 
	adPack, deriv=TRUE, 
	package = 'adlaplace') {

	Nbeta = nrow(data$XTp)
	config_inner = config
	config_inner$beta = x[seq(1, len=Nbeta)]
	config_inner$theta = x[seq(Nbeta+1, len=length(x) - Nbeta)]
	Sgamma1 = seq(Nbeta+1, len=length(start_gamma))

	if(length(start_gamma) != nrow(data$ATp)) {
		warning("start_gamma is the wrong size")
	}

	if((1+length(unique(data$map))) != length(config_inner$theta)) {
		warning("x is the wrong size")
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

	if(FALSE) {
		# check hessian chol

	checkHessian1 = 
		inner_res$cholHessian$L %*% Matrix::Diagonal(length(inner_res$cholHessian$D), inner_res$cholHessian$D) %*%
			Matrix::t(inner_res$cholHessian$L) 
		checkHessian2 = checkHessian1[1+inner_res$cholHessian$P, inner_res$cholHessian$P+1]
		check3 = (inner_res$hessian - checkHessian2)
	quantile(check3@x)

	}

	inner_res$parameters = x
	inner_res$fullParameters = c(config_inner$beta, inner_res$solution, config_inner$theta)

	if(!deriv & is.null(adPack)) {
		return(inner_res)
	}

	inner_res$cholHessian$halfH = reformatChol(inner_res$cholHessian)
	inner_res$cholHessian$Hinv = Matrix::crossprod(inner_res$cholHessian$halfH) 

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

	whichColumnsByGroup1 = lapply(
		config$group_inner, function(xx, refmat) {
			linvHere = refmat[1+xx$grad, ,drop=FALSE]
			which(diff(linvHere@p)>0)-1L
		}, 
		refmat = inner_res$cholHessian$halfH
	)

	whichColumnsByGroup = Matrix::sparseMatrix(
		i = unlist(whichColumnsByGroup1),
		j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
		index1=FALSE,
		dims = c(ncol(inner_res$cholHessian$halfH), length(whichColumnsByGroup1))
	)

	theTrace = getExportedValue(package, "traceHinvT")(
		inner_res$fullParameters, 
		inner_res$cholHessian$halfH, 
		whichColumnsByGroup,
		config,
		adPack
	)


	dU = -result$inner$cholHessian$Hinv %*% result$outer$hessian[Sgamma1, -Sgamma1]

	result$deriv = data.frame(
		dDetUpart = as.vector(theTrace[Sgamma1] %*% dU),
		dDetTpart = theTrace[-Sgamma1])
	result$deriv$gradTheta = result$outer$grad[-Sgamma1]  
	result$deriv$gradU = as.vector(result$outer$grad[Sgamma1] %*% dU)
	result$deriv$dDet = result$deriv$dDetUpart + result$deriv$dDetTpart
	result$deriv$dL = result$deriv$dDet + result$deriv$gradU + result$deriv$gradTheta

	result$dLogLik = result$deriv$dL

return(result)
}

