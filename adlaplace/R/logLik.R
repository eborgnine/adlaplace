#' Log-likelihood with inner Laplace optimization
#'
#' Evaluates the (profiled) log-likelihood for a hierarchical model by solving an
#' inner optimization problem for the latent vector \code{gamma} (e.g. random
#' effects) given outer parameters \code{x} (typically \code{beta} and \code{theta}).
#' Optionally computes derivatives via automatic differentiation using an AD
#' “tape” object (\code{adPack}).
#'
#' The function delegates the inner optimization to \code{inner_opt()} from the
#' selected backend package (by default \pkg{adlaplace}), and (when \code{deriv=TRUE})
#' computes derivatives.
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
#'   parameter \code{gamma}. Defaults to \code{config$start_gamma}.
#' @param control List of control parameters passed to the backend inner optimizer.
#'   (e.g. \code{report.level}, \code{report.freq}).
#' @param adPack Optional AD object returned by the backend \code{getAdFun()}.
#'   If missing and \code{deriv=TRUE}, it will be constructed automatically.
#' @param deriv Logical; if \code{TRUE}, compute derivatives.
#' @param package Character scalar naming the backend package to use for
#'   \code{getAdFun()} and \code{inner_opt()}. Defaults to the first element of
#'   \code{c(config$package, "adlaplace")}.
#'
#' @details
#' The parameter vector \code{x} is split into \code{beta} and \code{theta} and
#' inserted into \code{config} (as \code{config_inner$beta} and
#' \code{config_inner$theta}) before calling the backend inner optimizer.
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
#'   \item{deriv}{(Only if \code{deriv=TRUE}) derivative output.}
#'   \item{dLogLik}{(Only if \code{deriv=TRUE}) shortcut to \code{result$deriv$deriv$dL}.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # x <- c(beta, theta)
#' # out <- logLik(x, data, config, deriv = TRUE)
#' # out$minusLogLik
#' # out$dLogLik
#' }
#'

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

	if(any(config$verbose)) {
		cat("logLik using package ", package, "for objective funcion\n")
	}

	if(missing(adPack)) {
		if(deriv) {
			adPack = getExportedValue(package, "getAdFun")(data, config, inner=FALSE)
		} else {
			adPack = NULL
		}
	}


	inner_res = try(getExportedValue(package, "inner_opt")(
		start_gamma,
		data=data, 
		config=config_inner, 
		control=control,
		adPackFull = adPack))
	if(any(class(inner_res) == 'try-error')) {
		cat("resetting starting values to all zero\n")
		cat("theta ", paste(x, collapse=" "), "\n")
		inner_res = try(getExportedValue(package, "inner_opt")(
		rep(0.0, length(start_gamma)),
		data=data, 
		config=config_inner, 
		control=control,
		adPackFull = adPack))		
	}

	if(any(config$verbose)) {
		cat("done inner opt\n")
	}

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

