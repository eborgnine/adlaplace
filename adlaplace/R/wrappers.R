#' Outer objective and gradient wrappers
#'
#' Convenience wrappers around \code{\link{log_lik_laplace}} for use with outer optimizers
#' that expect \code{fn} / \code{gr} callbacks. Both functions solve (or warm-start)
#' the inner problem over \code{gamma} via \code{\link{log_lik_laplace}} and update a mutable
#' \code{cache} environment with the latest inner solution.
#'
#' \describe{
#' \item{\code{outer_fn()}}{Returns the scalar objective negative log likelihood.}
#' \item{\code{outer_gr()}}{Returns the gradient.}
#' }
#'
#' @param x Numeric outer parameter vector \code{c(beta, theta)}.
#' @param config Configuration list passed to \code{\link{log_lik_laplace}}.
#' @param ad_fun \code{ad_fun} object from \code{\link{ad_fun}}.
#' @param ... Additional arguments forwarded to \code{\link{log_lik_laplace}}
#'   (for example \code{data} or \code{package}).
#' @param control_inner A list of control options forwarded to the \code{control}
#'   argument of \code{\link{log_lik_laplace}} for the inner optimization.
#' @param cache An \code{\link[base]{environment}} containing starting values for the inner
#'   optimization. It should contain \code{gamma}. Both functions update
#'   \code{cache$gamma} to the latest \code{gamma} solution.
#'
#' @return
#' \itemize{
#' \item \code{outer_fn}: a numeric scalar \code{neg_log_lik}.
#' \item \code{outer_gr}: gradient of \code{neg_log_lik} (same sign as \code{outer_fn}).
#' }
#'
#' @seealso \code{\link{log_lik_laplace}}
#'
#' @examples
#' \dontrun{
#' cache <- new.env(parent = emptyenv())
#' cache$gamma <- rep(0, nrow(data$ATp))
#' ad_fun <- ad_fun(data, config)
#'
#' val <- outer_fn(x = x0, data = data, config = config, cache = cache, ad_fun = ad_fun)
#' gr <- outer_gr(x = x0, data = data, config = config, cache = cache, ad_fun = ad_fun)
#' }
#'
#' @name outer_optim_wrappers
#' @rdname outer_optim_wrappers
#' @export
outer_fn <- function(x, config, cache, ad_fun, control_inner = list(), ...) {
    assign("last_par_fn", x, cache)
    if (is.null(config$gamma)) {
        stop("outer_fn requires config$gamma")
    }
    if (is.null(cache$gamma) || length(cache$gamma) != length(config$gamma)) {
        cache$gamma <- config$gamma
    }

    result <- adlaplace::inner_opt(
        parameters = x,
        gamma = cache$gamma,
        config = config,
        ad_fun = ad_fun,
        control = control_inner,
        deriv = FALSE
    )

    assign("gamma", result$opt$solution, cache)
    result$neg_log_lik
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr <- function(x, config, cache, ad_fun, control_inner = list(), ...) {
    assign("last_par_gr", x, cache)
    if (is.null(config$gamma)) {
        stop("outer_gr requires config$gamma")
    }
    if (is.null(cache$gamma) || length(cache$gamma) != length(config$gamma)) {
        cache$gamma <- config$gamma
    }

    config_inner <- config

    result <- adlaplace::log_lik_laplace(
        x = x, config = config_inner,
        gamma = cache$gamma,
        control = control_inner,
        ad_fun = ad_fun,
        deriv = TRUE, ...
    )
    assign("gamma", result$opt$solution, cache)
    result$grad
}
