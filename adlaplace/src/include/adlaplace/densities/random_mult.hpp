#ifndef ADLAPLACE_RANDOM_MULT_HPP
#define ADLAPLACE_RANDOM_MULT_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/creators/ad_data.hpp"

CppAD::vector<CppAD::AD<double>> random_mult(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  SEXP precision,
  const Rcpp::List& config);

#endif
