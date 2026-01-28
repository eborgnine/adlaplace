#' Outer objective and gradient wrappers
#'
#' Convenience wrappers around \code{\link{logLik}} for use with outer optimizers
#' that expect \code{fn} / \code{gr} callbacks. Both functions solve (or warm-start)
#' the inner problem over \code{gamma} via \code{\link{logLik}} and update a mutable
#' \code{cache} environment with the latest inner solution.
#'
#' \describe{
#' \item{\code{outer_fn()}}{Returns the scalar objective negative log likelihood.}
#' \item{\code{outer_gr()}}{Returns the gradient.}
#' }
#'
#' @param ... Arguments forwarded to \code{\link{logLik}} (e.g. \code{x}, \code{data},
#'   \code{config}, \code{adPack}, \code{package}, etc.).
#' @param control_inner A list of control options forwarded to the \code{control}
#'   argument of \code{\link{logLik}} for the inner optimization.
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
#' @seealso \code{\link{logLik}}
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
	result = adlaplace::logLik(..., control = control_inner, start_gamma = cache$start_gamma, deriv=FALSE)
	assign('start_gamma', result$inner$solution,  cache)
	result$minusLogLik
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr = function(..., control_inner = list(), cache) {
	result = adlaplace::logLik(..., control = control_inner, start_gamma = cache$start_gamma, deriv=TRUE)
	assign('start_gamma', result$inner$solution, cache)
	result$dLogLik
}