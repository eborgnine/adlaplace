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
#' @param cache An \code{\link[base]{environment}} used to warm-start and store the inner
#'   \code{gamma} solution. If \code{cache$gamma} is missing or the wrong length,
#'   it is initialized from \code{config$gamma} (if present) or zeros of length
#'   \code{ad_fun@sizes["gamma"]}. Both functions update \code{cache$gamma} after
#'   each evaluation.
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
outer_fn <- function(
  x, config, cache, ad_fun, control_inner = list(), ...
) {
    assign("last_par_fn", x, cache)
    num_gamma <- as.integer(ad_fun@sizes["gamma"])
    cache$gamma <- resolve_gamma_start(config, cache, num_gamma)

    result <- adlaplace::inner_opt(
        parameters = x,
        gamma = cache$gamma,
        ad_fun = ad_fun,
        control = control_inner,
        deriv = FALSE,
        verbose = isTRUE(config$verbose)
    )

    assign("gamma", result$opt$solution, cache)
    result$neg_log_lik
}

#' @rdname outer_optim_wrappers
#' @export
outer_gr <- function(
  x, config, cache, ad_fun, control_inner = list(), ...
) {
    assign("last_par_gr", x, cache)
    num_gamma <- as.integer(ad_fun@sizes["gamma"])
    cache$gamma <- resolve_gamma_start(config, cache, num_gamma)

    result <- adlaplace::log_lik_laplace(
        x = x, config = config,
        gamma = cache$gamma,
        control = control_inner,
        ad_fun = ad_fun,
        deriv = TRUE, ...
    )
    assign("gamma", result$opt$solution, cache)
    result$grad
}
