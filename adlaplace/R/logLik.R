#' @export
logLik = function(x, data, config, 
	start_gamma = config$start_gamma, 	
	control = list(report.level=4, report.freq=1), 
	adPack, deriv=TRUE) {

	Nbeta = nrow(data$XTp)
	config_inner = config
	config_inner$beta = x[seq(1, len=Nbeta)]
	config_inner$theta = x[seq(Nbeta+1, len=length(x) - Nbeta)]
	Sgamma1 = seq(Nbeta+1, len=length(start_gamma))


	inner_res = adlaplace::inner_opt(
		start_gamma,
		data=data, config=config_inner, 
		control=control)


	inner_res$parameters = x
	inner_res$fullParameters = c(config_inner$beta, inner_res$solution, config_inner$theta)

	if(!deriv) {
		return(inner_res)
	}

	if(missing(adPack)) {
		adPack = adlaplace::getAdFun(data, config, inner=FALSE)
	}

	result = c(
		list(
			inner = inner_res[grep("[pP]arameters", names(inner_res), invert=TRUE)],
			outer = list(hessian = adlaplace::hessian(inner_res$fullParameters, adPack, config),
				grad = adlaplace::grad(inner_res$fullParameters, adPack, config)
			)
		),
		inner_res[grep("[pP]arameters|minusLogLik", names(inner_res))]
	)


	Linv = Matrix::solve(inner_res$cholHessian$L)
	#Linv[theseq, theseq]
		# H = Pt L D Lt P
		# Hinv = Pt LinvT Dinv Linv P
	halfDinv = Matrix::Diagonal(length(inner_res$cholHessian$D), (inner_res$cholHessian$D)^(-0.5))
	halfH = (halfDinv %*% Linv)[,1+inner_res$cholHessian$P] 
	LinvPt = Matrix::t(halfH)
	Hinv =  halfH %*% LinvPt

	linvL = as(LinvPt, 'lMatrix')

	whichColumnsByGroup1 = #parallel::mc
	lapply(
		config$group_inner, function(xx) {
			linvHere = linvL[1+xx$grad, ,drop=FALSE]
			which(diff(linvHere@p)>0)-1L
		}#, mc.cores = max(c(config$num_threads, 1), na.rm=TRUE)
	)
	whichColumnsByGroup = Matrix::sparseMatrix(
		i = unlist(whichColumnsByGroup1),
		j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
		index1=FALSE,
		dims = c(nrow(Hinv), length(whichColumnsByGroup1))
	)

	theTrace = adlaplace::traceHinvT(
		inner_res$fullParameters, 
		LinvPt, 
		whichColumnsByGroup,
		config,
		adPack
	)


	dU = -Hinv %*% result$outer$hessian[Sgamma1, -Sgamma1]

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

