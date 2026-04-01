#' @useDynLib hpolcc, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom adlaplace design precision theta_info beta_info random_info
#' @importFrom adlaplace f
#' 
.my_beta_init <- 0
.my_beta_lower <- -Inf
.my_beta_upper <- Inf
.my_beta_parscale <- 1

.my_theta_init <- 0.02
.my_theta_lower <- 1e-9
.my_theta_upper <- Inf
.my_theta_parscale <- 1

