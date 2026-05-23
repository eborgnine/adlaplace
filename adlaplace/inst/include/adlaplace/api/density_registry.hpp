#ifndef ADLAPLACE_DENSITY_REGISTRY_HPP
#define ADLAPLACE_DENSITY_REGISTRY_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <string>

#include "adlaplace/creators/ad_model.hpp"

using LogDensObsFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_model& model,
  const Rcpp::List& config,
  size_t Dgroup);

using LogDensSingleDataFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_model& model,
  const Rcpp::List& config);

using LogDensSingleRandomDiagFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_model& model,
  const Rcpp::List& precision,
  const Rcpp::List& config);

void register_log_dens_obs(const std::string& name, LogDensObsFn fn);
void register_log_dens_single_data(const std::string& name, LogDensSingleDataFn fn);
void register_log_dens_single_random_diag(const std::string& name, LogDensSingleRandomDiagFn fn);

LogDensObsFn resolve_log_dens_obs(const std::string& name);
LogDensSingleDataFn resolve_log_dens_single_data(const std::string& name);
LogDensSingleRandomDiagFn resolve_log_dens_single_random_diag(const std::string& name);

bool log_dens_single_uses_random_diag(const std::string& name);

void register_adlaplace_default_densities();

#endif  // ADLAPLACE_DENSITY_REGISTRY_HPP
