#' @useDynLib adlaplace, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @import methods
#' @import RCppAD
NULL

loadNamespace("Matrix")

.onLoad <- function(libname, pkgname) {
  if (isTRUE(has_openmp())) {
    warm_openmp_runtime()
  }
  invisible(NULL)
}

#' Raw AD handle (external pointer)
#'
#' S4 class for a C++ AD backend handle. Objects are created by
#' \code{\link{ad_pack_ptr}}, combined with \code{\link{c.ad_pack_ptr}}, and
#' typically passed to \code{\link{ad_pack}} to attach Hessian templates.
#'
#' @name ad_pack_ptr-class
#' @aliases ad_pack_ptr-class
#' @exportClass ad_pack_ptr
#' @seealso \code{\link{ad_pack_ptr}}, \code{\link{ad_pack}},
#'   \code{\link{clone_ad_pack_ptr}}
setClass("ad_pack_ptr", contains = "externalptr")

#' @keywords internal
setOldClass("ad_pack_ptr")

.my_beta_init <- 0
.my_beta_lower <- -Inf
.my_beta_upper <- Inf
.my_beta_parscale <- 1

.my_theta_init <- 0.02
.my_theta_lower <- 1e-9
.my_theta_upper <- Inf
.my_theta_parscale <- 1

#' @keywords internal
ad_shard_label <- function(shard, name = NULL) {
  if (!is.null(name) && length(name) == 1L && nzchar(name)) {
    return(name)
  }
  paste0(shard@ad_kind, "/", shard@density)
}

#' @keywords internal
init_from_info_block <- function(block) {
  if (is.null(block) || !is.data.frame(block) || nrow(block) == 0L) {
    return(numeric(0))
  }
  as.numeric(block$init)
}

#' @keywords internal
is_model_data_bundle <- function(x) {
  is.list(x) &&
    all(c("term_data", "observations", "random", "parameters") %in% names(x)) &&
    is.list(x$observations) &&
    is.list(x$random) &&
    is.list(x$parameters)
}
