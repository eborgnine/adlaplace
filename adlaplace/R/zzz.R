#' @useDynLib adlaplace, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @import methods
NULL

loadNamespace("Matrix")

#' @exportClass ad_fun_ptr
setClass("ad_fun_ptr", contains = "externalptr")

#' @keywords internal
setOldClass("ad_fun_ptr")

.my_beta_init <- 0
.my_beta_lower <- -Inf
.my_beta_upper <- Inf
.my_beta_parscale <- 1

.my_theta_init <- 0.02
.my_theta_lower <- 1e-9
.my_theta_upper <- Inf
.my_theta_parscale <- 1

#' @keywords internal
ad_fun_shard_label <- function(shard, name = NULL) {
  if (!is.null(name) && length(name) == 1L && nzchar(name)) {
    return(name)
  }
  paste0(shard@ad_kind, "/", shard@ad_fun)
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
    all(c("data", "observations", "random", "parameters") %in% names(x)) &&
    is.list(x$observations) &&
    is.list(x$random) &&
    is.list(x$parameters)
}
