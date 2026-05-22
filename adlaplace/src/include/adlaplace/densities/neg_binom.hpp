#ifndef ADLAPLACE_NEG_BINOM_HPP
#define ADLAPLACE_NEG_BINOM_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/creators/rviews.hpp"

CppAD::vector<CppAD::AD<double>> neg_binom_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const Rcpp::List& data,
  const Rcpp::List& config,
  size_t Dgroup);

CppAD::vector<CppAD::AD<double>> neg_binom_extra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const Rcpp::List& data,
  const Rcpp::List& config);

#endif
