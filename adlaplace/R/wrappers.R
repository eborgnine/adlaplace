#' @export
outer_fn = function(..., control_inner = list(), cache) {
	result = adlaplace::logLik(..., control = control_inner, start_gamma = cache$start_gamma, deriv=FALSE)
	assign('start_gamma', result$solution,  cache)
	assign('x', result$parameters,  cache)
	result$minusLogLik
}

#' @export
outer_gr = function(..., control_inner = list(), cache) {
	result = adlaplace::logLik(..., control = control_inner, start_gamma = cache$start_gamma, deriv=TRUE)
	assign('start_gamma', result$inner$solution, cache)
	result$dLogLik
}