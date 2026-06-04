#' @useDynLib adlaplaceExample, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.onLoad <- function(libname, pkgname) {
  register_example_densities()
}
