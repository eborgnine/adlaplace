#ifndef ADLAPLACE_NBINOM_HPP
#define ADLAPLACE_NBINOM_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/creators/ad_data.hpp"

CppAD::vector<CppAD::AD<double>> nbinom_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config,
  size_t Dgroup);

CppAD::vector<CppAD::AD<double>> nbinom_extra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config);

#endif
