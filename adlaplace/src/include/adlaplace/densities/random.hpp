#ifndef ADLAPLACE_RANDOM_HPP
#define ADLAPLACE_RANDOM_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

CppAD::vector<CppAD::AD<double>> random_diagonal(
  const CppAD::vector<CppAD::AD<double>>& x,
  const Rcpp::List& data,
  const Rcpp::List& config);

#endif
