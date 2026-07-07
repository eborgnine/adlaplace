#' @describeIn adlaplace_cpp Accept \code{ad_fun} S4 (uses \code{@ptr}).
#' @export
joint_log_dens <- function(ad_fun, x, shards = NULL, negative = TRUE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call(`_adlaplace_joint_log_dens`, ptr, x, shards, negative)
}

#' @describeIn adlaplace_cpp Accept \code{ad_fun} S4 (uses \code{@ptr}).
#' @export
grad <- function(ad_fun, x, shards = NULL, inner = FALSE, negative = TRUE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call(`_adlaplace_grad`, ptr, x, shards, inner, negative)
}

#' @describeIn adlaplace_cpp Accept \code{ad_fun} S4 (uses \code{@ptr}).
#' @export
hessian <- function(ad_fun, x, shards = NULL, inner = FALSE, verbose = FALSE, negative = TRUE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call(`_adlaplace_hessian`, ptr, x, shards, inner, verbose, negative)
}

#' @describeIn adlaplace_cpp Accept \code{ad_fun} S4 (uses \code{@ptr}).
#' @param verbose Logical; if \code{TRUE}, print threads, shards, and CppAD session phases.
#' @export
trace_hinv_t <- function(ad_fun, x, LinvPt, LinvPtColumns, verbose = FALSE) {
  ptr <- if (isS4(ad_fun) && methods::.hasSlot(ad_fun, "ptr")) ad_fun@ptr else ad_fun
  .Call("_adlaplace_trace_hinv_t", ptr, x, LinvPt, LinvPtColumns, verbose)
}

#' @describeIn adlaplace_cpp Accept \code{ad_fun} S4 (requires Hessian templates).
#' @param parameters Numeric vector of length \code{Nbeta + Ntheta}.
#' @param gamma Numeric vector of length \code{Ngamma}.
#' @export
fun_obj_fdfh <- function(ad_fun, parameters, gamma, inner = TRUE, verbose = FALSE) {
  .Call(`_adlaplace_fun_obj_fdfh`, parameters, gamma, ad_fun, inner, verbose)
}
