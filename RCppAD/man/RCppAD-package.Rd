\name{RCppAD-package}
\alias{RCppAD}
\alias{RCppAD-package}
\docType{package}
\title{CppAD C++ Header Files for Automatic Differentiation}
\description{
  Provides the CppAD C++ header library for automatic differentiation,
  for use by R packages via \code{LinkingTo: RCppAD}.
}
\details{
  After depending on this package, include CppAD from C++ with
  \code{#include <cppad/cppad.hpp>}.

  Headers are vendored with optional backends (ColPack, Eigen ADOLC,
  IPOPT) disabled. The error handler uses \code{REprintf} /
  \code{Rf_error} instead of \code{std::cerr} / \code{std::exit}.
}
\author{EPL-2.0 | GPL (>= 2)}
\keyword{package}
