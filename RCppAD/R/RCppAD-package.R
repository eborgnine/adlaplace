#' RCppAD: CppAD C++ headers for R packages
#'
#' This package ships the CppAD header library under
#' \code{inst/include/cppad/}. Packages that need CppAD should add
#' \code{LinkingTo: RCppAD} (and usually \code{Imports: RCppAD}) in
#' \code{DESCRIPTION}, then \code{#include <cppad/cppad.hpp>} from C++.
#'
#' Refresh upstream headers with \code{tools/update-cppad.sh} (not run at
#' install time).
#'
#' @name RCppAD-package
#' @keywords internal
"_PACKAGE"
