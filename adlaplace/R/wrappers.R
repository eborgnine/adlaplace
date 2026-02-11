#' Outer objective and gradient wrappers
#'
#' Convenience wrappers around \code{\link{logLikLaplace}} for use with outer optimizers
#' that expect \code{fn} / \code{gr} callbacks. Both functions solve (or warm-start)
#' the inner problem over \code{gamma} via \code{\link{logLikLaplace}} and update a mutable
#' \code{cache} environment with the latest inner solution.
#'
#' \describe{
#' \item{\code{outer_fn()}}{Returns the scalar objective negative log likelihood.}
#' \item{\code{outer_gr()}}{Returns the gradient.}
#' }
#'
#' @param ... Arguments forwarded to \code{\link{logLikLaplace}} (e.g. \code{x}, \code{data},
#'   \code{config}, \code{adFun}, \code{package}, etc.).
#' @param control_inner A list of control options forwarded to the \code{control}
#'   argument of \code{\link{logLikLaplace}} for the inner optimization.
#' @param cache An \code{\link[base]{environment}} containing starting values for the inner
#'   optimization. It should contain \code{start_gamma}. Both functions update
#'   \code{cache$start_gamma} to the latest \code{gamma} solution.
#'
#' @return
#' \itemize{
#' \item \code{outer_fn}: a numeric scalar (objective value).
#' \item \code{outer_gr}: a numeric vector (gradient w.r.t. outer parameters).
#' }
#'
#' @seealso \code{\link{logLikLaplace}}
#'
#' @examples
#' \dontrun{
#' cache <- new.env(parent = emptyenv())
#' cache$start_gamma <- rep(0, nrow(data$ATp))
#'
#' val <- outer_fn(x = x0, data = data, config = config, cache = cache)
#' gr  <- outer_gr(x = x0, data = data, config = config, cache = cache)
#' }
#'
#' @name outer_optim_wrappers
#' @rdname outer_optim_wrappers
#' @export
outer_fn = function(x, config, cache, adFun, control_inner = list(), ...) {
	if(is.null(cache$start_gamma)) {
		warning("starting values cache$start_gamma not provided")
	}
	result = adlaplace::logLikLaplace(
		x=x, config=config,
		start_gamma = cache$start_gamma,
		control = control_inner,
		adFun = adFun, 
		deriv = FALSE, ...
	)
	assign('start_gamma', result$opt$solution, cache)
	result$fval
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr = function(x, config, cache, adFun, control_inner = list(), ...) {
	result = adlaplace::logLikLaplace(
		x=x, config=config,
		start_gamma = cache$start_gamma,
		control = control_inner,
		adFun = adFun, 
		deriv = TRUE, ...
	)
	assign('start_gamma', result$opt$solution, cache)
	result$grad
}
