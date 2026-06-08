#ifndef ADLAPLACE_LOG_DENS_FN_HPP
#define ADLAPLACE_LOG_DENS_FN_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <cstddef>

#include "adlaplace/creators/ad_data.hpp"

using LogDensObsFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config,
  size_t Dgroup);

using LogDensSingleDataFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config);

using LogDensSingleRandomDiagFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  SEXP precision,
  const Rcpp::List& config);

#endif  // ADLAPLACE_LOG_DENS_FN_HPP
