#ifndef ADFUN_SINGLE_HPP
#define ADFUN_SINGLE_HPP

#include "adlaplace/api/density_registry.hpp"
#include "adlaplace/creators/adfun_common.hpp"
#include "adlaplace/creators/rviews.hpp"

inline GroupPack build_ad_fun_single(
  const Rcpp::List& data,
  const Rcpp::List& config,
  const std::string& single_name) {

  const Config cfg(config);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_single " << single_name << "\n";
  }

  CPPAD_TESTVECTOR(double) ad_params_G(cfg.Nparams);
  for (size_t d = 0; d < cfg.Nparams; ++d) {
    ad_params_G[d] = cfg.params[d];
  }

  CppAD::vector<CppAD::AD<double>> ad_params(cfg.Nparams);
  for (size_t d = 0; d < cfg.Nparams; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  CppAD::Independent(ad_params);

  CppAD::vector<CppAD::AD<double>> result_here;
  if (log_dens_single_uses_random_diag(single_name)) {
    LogDensSingleRandomDiagFn log_dens = resolve_log_dens_single_random_diag(single_name);
    result_here = log_dens(ad_params, data, config);
  } else {
    LogDensSingleDataFn log_dens = resolve_log_dens_single_data(single_name);
    result_here = log_dens(ad_params, data, config);
  }

  CppAD::ADFun<double> fun(ad_params, result_here);

  GroupPack pack;
  pack.fun = std::move(fun);
  adpack_sparsity(ad_params_G, cfg.Sgamma, pack, cfg.verbose);
  return pack;
}

#endif
