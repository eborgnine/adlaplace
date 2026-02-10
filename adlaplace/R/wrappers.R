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
outer_fn = function(..., control_inner = list(), cache) {
	result = adlaplace::logLikLaplace(..., control = control_inner, start_gamma = cache$start_gamma)
	assign('start_gamma', result$inner$solution,  cache)
	result$minusLogLik
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr = function(..., control_inner = list(), cache) {
	args <- list(...)
	if(is.null(args$x)) {
		stop("outer_gr requires argument 'x' in ...")
	}
	x0 <- args$x
	eps <- 1e-6

	base <- do.call(
		adlaplace::logLikLaplace,
		c(args, list(control = control_inner, start_gamma = cache$start_gamma))
	)
	assign('start_gamma', base$inner$solution, cache)

	gr <- numeric(length(x0))
	for(i in seq_along(x0)) {
		xp <- x0; xm <- x0
		xp[i] <- xp[i] + eps
		xm[i] <- xm[i] - eps

		args_p <- args; args_m <- args
		args_p$x <- xp
		args_m$x <- xm

		fp <- do.call(
			adlaplace::logLikLaplace,
			c(args_p, list(control = control_inner, start_gamma = base$inner$solution))
		)$minusLogLik
		fm <- do.call(
			adlaplace::logLikLaplace,
			c(args_m, list(control = control_inner, start_gamma = base$inner$solution))
		)$minusLogLik
		gr[i] <- (fp - fm)/(2 * eps)
	}
	gr
}
