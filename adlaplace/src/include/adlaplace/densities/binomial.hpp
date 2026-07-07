#ifndef ADLAPLACE_BINOMIAL_HPP
#define ADLAPLACE_BINOMIAL_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/creators/ad_data.hpp"

CppAD::vector<CppAD::AD<double>> binomial_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config,
  size_t Dgroup);

#endif
