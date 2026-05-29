#' @rdname innerOpt
#' @export
inner_opt <- function(parameters, gamma, ad_fun, control = list(), deriv = FALSE, verbose = FALSE) {
  .Call(`_adlaplace_inner_opt`, parameters, gamma, ad_fun, control, deriv, verbose)
}
