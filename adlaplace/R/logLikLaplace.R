#' Log-likelihood with inner Laplace optimization
#'
#' Evaluates the (profiled) log-likelihood for a hierarchical model by solving an
#' inner optimization problem for the latent vector \code{gamma} (e.g. random
#' effects) given outer parameters \code{x} (typically \code{beta} and \code{theta}).
#' The function delegates the inner optimization to \code{inner_opt()} from the
#' selected backend package (by default \pkg{adlaplace}).
#'
#' @param x Numeric vector of outer parameters. The first \code{Nbeta} entries are
#'   interpreted as \code{beta}; the remainder are interpreted as \code{theta}
#'   (with \code{Nbeta = nrow(data$XTp)}).
#' @param data A list containing model matrices and metadata used by the backend.
#'   At minimum this function assumes components \code{XTp}, \code{ATp}, and
#'   \code{map} exist.
#' @param config A list of configuration options passed to the backend. This
#'   function uses \code{config$package} and \code{config$verbose} if present.
#' @param start_gamma Optional numeric vector of starting values for the inner
#'   parameter \code{gamma}. Defaults to \code{config$gamma}.
#' @param control List of control parameters passed to the backend inner optimizer.
#'   (e.g. \code{report.level}, \code{report.freq}).
#' @param adPack Optional AD object returned by the backend \code{getAdFun()}.
#'   This is a single backend handle (no separate inner/outer handles). If
#'   missing, it will be constructed automatically.
#' @param package Character scalar naming the backend package to use for
#'   \code{getAdFun()} and \code{inner_opt()}. Defaults to the first element of
#'   \code{c(config$package, "adlaplace")}.
#'
#' @details
#' The parameter vector \code{x} is split into \code{beta} and \code{theta} and
#' inserted into \code{config} (as \code{config_inner$beta} and
#' \code{config_inner$theta}) before calling the backend inner optimizer.
#'
#' Current backends use a single AD handle. This function passes that handle to
#' \code{inner_opt(..., adPack = adPack)}.
#'
#' The returned list contains both “inner” outputs (from the inner optimization)
#' and “outer” outputs (quantities associated with the outer objective), as well
#' as convenience fields such as \code{full_parameters}.
#'
#' @return A list with components including:
#' \describe{
#'   \item{inner}{List of inner-optimization outputs (excluding parameter bookkeeping).}
#'   \item{outer}{List of outer quantities (backend-dependent).}
#'   \item{minusLogLik}{Numeric scalar (backend-dependent sign convention).}
#'   \item{parameters}{The input outer parameter vector \code{x}.}
#'   \item{full_parameters}{Concatenation of \code{beta}, optimized \code{gamma},
#'     and \code{theta} (names are set when possible).}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # x <- c(beta, theta)
#' # out <- logLikLaplace(x, data, config)
#' # out$minusLogLik
#' }
#'

#' @export
logLikLaplace = function(
	x, config, 
	start_gamma, 	
	control = list(report.level=4, report.freq=1), 
	adPack, data, 
	package = c(config$package, 'adlaplace')[1],
	deriv = FALSE
) {

	Nbeta = length(config$beta)
	Ntheta = length(config$theta)
	config_inner = config
	config_inner$beta = x[seq.int(1, len=Nbeta)]
	config_inner$theta = x[seq.int(Nbeta+1, len=Ntheta)]
	if(!missing(start_gamma)) {
		config_inner$gamma = start_gamma
		if(length(config$gamma) != length(config_inner$gamma)) {
			warning("gamma is the wrong length")
		}
	} 

	if(missing(adPack)) {
		if(missing(data)) {
			stop("at least one of data and adPack must be supplied")
		}
		adPack = getExportedValue(package, "getAdFun")(data, config)
	}


	if(any(config$verbose)) {
		cat("logLikLaplace using package ", package, "for objective funcion\n")
	}

	inner_res = try(adlaplace::inner_opt(
		x, 
		config_inner$gamma,
		config=config_inner, 
		control=control,
		adPack = adPack))

	if(any(class(inner_res) == 'try-error')) {
		cat("resetting starting values to all zero\n")
		cat("theta ", paste(x, collapse=" "), "\n")
		config_inner$gamma = rep(0.0, length(config_inner$gamma))
		inner_res = try(inner_opt(
			config=config_inner, 
			control=control,
			adPack = adPack))
	}

	if(any(config$verbose)) {
		cat("done inner opt\n")
	}

	if(!missing(data)) {
		names(inner_res$solution) = rownames(data$ATp)
	} else {
		names(inner_res$solution) = names(config$gamma)
	}


	result = list(
		minusLogLik = inner_res$minusLogLik,
		parameters = x,
		fullParameters =  c(config_inner$beta, inner_res$solution, config_inner$theta),
		hessian = list(
			H = do.call(Matrix::sparseMatrix, inner_res$hessian),
			L = do.call(Matrix::sparseMatrix, inner_res$cholHessian$L), 
			D = Matrix::Diagonal(length(inner_res$solution), inner_res$cholHessian$D),
			P = inner_res$cholHessian$P
		),
		opt = inner_res[grep("[hH]essian", names(inner_res), invert=TRUE)]
	)

	if(!missing(data)){
		try(names(result$fullParameters) <- c(
			rownames(data$XTp), 
			rownames(data$ATp), 
			colnames(data$map)
		))
	}




	if(!deriv) {
		return(result)
	}	
	result$deriv = logLikDeriv(
		fullParameters = result$fullParameters, 
		hessianPack = result$hessian,
		config, adPack)

	return(result)
}
