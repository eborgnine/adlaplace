#ifndef ADLAPLACE_RANDOM_DIAGONAL_HPP
#define ADLAPLACE_RANDOM_DIAGONAL_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "adlaplace/creators/ad_data.hpp"

CppAD::vector<CppAD::AD<double>> random_diagonal(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  SEXP precision,
  const Rcpp::List& config);

#endif
