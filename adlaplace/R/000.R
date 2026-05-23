#' @useDynLib adlaplace, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @import methods
NULL

loadNamespace("Matrix")

#' @export
.type_factor_levels <- c("fixed", "random", "response")

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
