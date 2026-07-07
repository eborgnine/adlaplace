#ifndef ADLAPLACE_GAUSSIAN_HPP
#define ADLAPLACE_GAUSSIAN_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/creators/ad_data.hpp"

CppAD::vector<CppAD::AD<double>> gaussian_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config,
  size_t Dgroup);

CppAD::vector<CppAD::AD<double>> gaussian_extra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config);

#endif
