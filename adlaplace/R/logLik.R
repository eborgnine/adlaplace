logLik = function(x, data, config, 
	start_gamma = config$start_gamma, 	
	control = list(report.level=4, report.freq=1), 
	adPack, deriv=TRUE) {

	Nbeta = nrow(data$XTp)
	config_inner = config
	config_inner$beta = x[seq(1, len=Nbeta)]
	config_inner$theta = x[seq(Nbeta+1, len=length(x) - Nbeta)]

	inner_res = adlaplace::inner_opt(
		start_gamma,
		data=data, config=config, 
		control=control)

	inner_res$parameters = x
	inner_res$fullParameters = c(config_inner$beta, inner_res$solution, config_inner$theta)

	if(!deriv) {
		return(inner_res)
	}

	# get full hessian
	# theseq = seq(to=length(inner_res$solution), len=8);inner_res$cholHessian$L[theseq, theseq]
	Linv = Matrix::solve(inner_res$cholHessian$L)
	#Linv[theseq, theseq]
		# H = Pt L D Lt P
		# Hinv = Pt LinvT Dinv Linv P
	halfDinv = Matrix::Diagonal(length(inner_res$cholHessian$D), (inner_res$cholHessian$D)^(-0.5))
	halfH = (halfDinv %*% Linv)[,1+inner_res$cholHessian$P] 
	LinvPt = Matrix::t(halfH)
	Hinv =  halfH %*% LinvPt
	# theseq= 1:5 
	# ( halfH %*% LinvPt)[theseq, theseq]
	#  Hinv[theseq, theseq]
	# Matrix::solve(inner_res$hessian)[theseq, theseq]
	linvL = as(LinvPt, 'lMatrix')
	# do this during setup?
	whichColumnsByGroup1 = parallel::mclapply(
		config$group, function(xx) {
			linvHere = linvL[1+xx$grad, ,drop=FALSE]
			which(diff(linvHere@p)>0)-1L
		}, mc.cores = max(c(config$num_threads, 1), na.rm=TRUE)
	)
	whichColumnsByGroup = Matrix::sparseMatrix(
		i = unlist(whichColumnsByGroup1),
		j = rep(seq(0, len=length(whichColumnsByGroup1)), unlist(lapply(whichColumnsByGroup1, length))),
		index1=FALSE,
		dims = c(nrow(Hinv), length(whichColumnsByGroup1))
	)

	result$fullHessian = hessian(result$fullParameters, data, config, adFunFull)
	result$fullGrad = grad(result$fullParameters, data, config, adFunFull)



		theTrace = hpolcc:::traceHinvT(
  result$fullParameters, 
  cholExpand$LinvPt, data, config, adFunFull
)

dU = -Hinv %*% result$fullHessian[Sgamma1, -Sgamma1]

		result$deriv = data.frame(
  dDetUpart = as.vector(theTrace[Sgamma1] %*% dU),
  dDetTpart = theTrace[-Sgamma1])
result$deriv$gradTheta = result$fullGrad[-Sgamma1]  
result$deriv$gradU = as.vector(result$fullGrad[Sgamma1] %*% dU)
result$deriv$dDet = result$deriv$dDetUpart + result$deriv$dDetTpart
result$deriv$dL = result$deriv$dDet + result$deriv$gradU + result$deriv$gradTheta

result$dLogLik = result$deriv$dL

	}
	return(inner_res)
}

