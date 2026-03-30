#' Log-likelihood with inner Laplace optimization
#'
#' Evaluates the (profiled) log-likelihood for a hierarchical model by solving an
#' inner optimization problem for the latent vector \code{gamma} (e.g. random
#' effects) given outer parameters \code{x} (typically \code{beta} and \code{theta}).
#' The function delegates the inner optimization to \code{inner_opt()} from the
#' selected backend package (by default \pkg{adlaplace}).
#'
#' @param x Numeric vector of outer parameters \code{c(beta, theta)} with length
#'   \code{length(config$beta) + length(config$theta)}.
#' @param config A list containing model dimensions/starting values and backend
#'   options. Must include \code{beta}, \code{gamma}, and \code{theta}; may also
#'   include \code{package}, \code{verbose}, and \code{num_threads}.
#' @param gamma Optional numeric vector of starting values for the inner
#'   parameter \code{gamma}. Defaults to \code{config$gamma}.
#' @param control List of control parameters passed to the backend inner optimizer.
#'   (e.g. \code{report.level}, \code{report.freq}).
#' @param adFun Optional AD object returned by the backend \code{getAdFun()}.
#'   This is a single backend handle (no separate inner/outer handles). If
#'   missing, it will be constructed automatically.
#' @param data Optional data list used to build \code{adFun} when \code{adFun}
#'   is not supplied.
#' @param package Character scalar naming the backend package to use for
#'   \code{getAdFun()} and \code{inner_opt()}. Defaults to the first element of
#'   \code{c(config$package, "adlaplace")}.
#' @param deriv Logical scalar. If \code{TRUE}, include derivative quantities in
#'   the output.
#'
#' @details
#' The parameter vector \code{x} is split into \code{beta} and \code{theta} and
#' inserted into \code{config} (as \code{config_inner$beta} and
#' \code{config_inner$theta}) before calling the backend inner optimizer.
#'
#' Current backends use a single AD handle. This function passes that handle to
#' \code{inner_opt(..., adFun = adFun)}.
#'
#' The inner objective is treated as negative joint log density; this function
#' returns both the Laplace-approximated log-likelihood (\code{logLik}) and its
#' negation (\code{fval}) for optimizer convenience.
#'
#' @return A list with components including:
#' \describe{
#'   \item{logLik}{Laplace-approximated log-likelihood at the optimized inner
#'   solution.}
#'   \item{fval}{\code{-logLik}.}
#'   \item{parameters}{The input outer parameter vector \code{x}.}
#'   \item{full_parameters}{Concatenation of \code{beta}, optimized \code{gamma},
#'   and \code{theta}.}
#'   \item{hessian}{List with \code{H} (outer Hessian sparse matrix) and
#'   \code{cholInner} (sparse LDL decomposition of inner Hessian).}
#'   \item{opt}{Inner optimizer outputs (excluding Hessian objects), including
#'   \code{solution}, \code{gradient}, \code{iterations}, \code{status}, and
#'   \code{trust.radius}.}
#'   \item{grad}{When \code{deriv=TRUE}: gradient currently returned with the
#'   package's outer-objective sign convention (negative log-likelihood).}
#'   \item{deriv, extra}{When \code{deriv=TRUE}: intermediate derivative pieces
#'   returned by \code{logLikDeriv()}.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # x <- c(beta, theta)
#' # out <- logLikLaplace(x = x, config = config, data = data, deriv = TRUE)
#' # out$logLik
#' # out$grad
#' }
#'

#' @export
logLikLaplace = function(
	x, config, 
	gamma, 	
	control = list(report.level=4, report.freq=1), 
	adFun, data, 
	package = c(config$package, 'adlaplace')[1],
	deriv = FALSE
) {

	Nbeta = length(config$beta)
	Ntheta = length(config$theta)
	Ngamma = length(config$gamma)
	Sgamma1 = seq.int(Nbeta+1, len=Ngamma)

	config_inner = config
	config_inner$beta = x[seq.int(1, length.out=Nbeta)]
	config_inner$theta = x[seq.int(Nbeta+1, length.out=Ntheta)]
	if(!missing(gamma)) {
		config_inner$gamma = gamma
		if(length(config$gamma) != length(config_inner$gamma)) {
			warning("gamma is the wrong length; resetting to config$gamma")
			config_inner$gamma = config$gamma
		}
	} 

	if(missing(adFun)) {
		if(missing(data)) {
			stop("at least one of data and adFun must be supplied")
		}
		adFun = adlaplace::getAdFun(data, config, package = package)
	} else {
		adfun_backend <- attr(adFun, "adlaplace.backend", exact = TRUE)
		if(!is.null(adfun_backend) && !identical(adfun_backend, package)) {
			stop(
				"adFun was built with backend package '", adfun_backend,
				"' but `package` is '", package, "'. ",
				"Rebuild with adlaplace::getAdFun(..., package = '", package, "')."
			)
		}
	}


	if(any(config$verbose)) {
		cat("logLikLaplace using package ", package, "for objective funcion\n")
	}

	Niter = 0;tryAgain=TRUE
	while(tryAgain & (Niter < 3) ) {
		Niter = Niter + 1
		inner_res = try(adlaplace::inner_opt(
			x, 
			config_inner$gamma,
			config=config_inner,
			control=control,
			adFun = adFun))

		tryAgain = any(class(inner_res) == 'try-error')
		if(!tryAgain) {
			tryAgain = sum(abs(inner_res$gradient[Sgamma1])) > 1
		}
		if(tryAgain) {
			cat("resetting starting values to all zero, ")
			cat("theta: ", paste(x, collapse=", "), "\n")
			config_inner$gamma = rep(0.0, length(config$gamma))
		}

	} # while
	if(any(class(inner_res) == 'try-error')) {
		stop("inner_opt failed in logLikLaplace: ", as.character(inner_res))
	}
	if(sum(abs(inner_res$gradient[Sgamma1])) > 1) {
		warning("inner_opt failed, large gradient")
	}
	if(any(config$verbose)) {
		cat("done inner opt\n")
	}

	Houter = do.call(Matrix::sparseMatrix, inner_res$hessian)
	Hinner = Houter[Sgamma1, Sgamma1]
	Hchol = Matrix::expand2(Matrix::Cholesky(Hinner, perm=TRUE, ldl=TRUE))
	
	halfLogDet = sum(log(Hchol$D@x))/2
	ONEHALFLOGTWOPI = 0.9189385332046727417803297364056176398613974736377834128171515404;

	logLik = -inner_res$fval - halfLogDet + Ngamma * ONEHALFLOGTWOPI;  


	result = list(
		logLik = logLik,
		fval = -logLik,
		parameters = x,
		full_parameters =  c(config_inner$beta, inner_res$solution, config_inner$theta),
		hessian = list(
			H = Houter,
			cholInner = Hchol
		),
		opt = inner_res[grep("[hH]essian", names(inner_res), invert=TRUE)]
	)

	if(!deriv) {
		return(result)
	}	

	theDeriv = logLikDeriv(
		full_parameters = result$full_parameters, 
		hessianPack = result$hessian,
		grad = inner_res$gradient,
		config, adFun)

	result$grad = -theDeriv$deriv$dL
	result$deriv = theDeriv$deriv
	result$extra = theDeriv$extra

	return(result)
}
