#ifndef ADLAPLACE_RANDOM_HPP
#define ADLAPLACE_RANDOM_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/creators/ad_model.hpp"

CppAD::vector<CppAD::AD<double>> random_diagonal(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_model& model,
  const Rcpp::List& precision,
  const Rcpp::List& config);

#endif
