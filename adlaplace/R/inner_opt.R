#' @rdname innerOpt
#' @export
inner_opt <- function(parameters, gamma, config, ad_fun, control = list(), deriv = FALSE) {
  .Call(`_adlaplace_inner_opt`, parameters, gamma, config, ad_fun, control, deriv)
}
