#' @importFrom methods setGeneric setMethod setOldClass as setClass new signature
#' @importFrom Matrix Matrix drop0 kronecker tcrossprod Diagonal sparseMatrix
#' @importFrom splines splineDesign
#' @importFrom adlaplace design precision theta_info beta_info random_info extra_ad_fun
#' @importFrom Rcpp evalCpp
#' @import RCppAD
#' @useDynLib adlaplaceFem, .registration = TRUE
NULL

# Allow S4 methods on data.frame (e.g. expand.grid output)
methods::setOldClass("data.frame")
