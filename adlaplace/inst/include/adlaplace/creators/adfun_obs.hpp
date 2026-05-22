#ifndef ADFUN_OBS_HPP
#define ADFUN_OBS_HPP

#include "adlaplace/api/density_registry.hpp"
#include "adlaplace/creators/adfun_common.hpp"
#include "adlaplace/creators/rviews.hpp"

inline size_t count_obs_groups(const Rcpp::List& data, const Rcpp::List& config) {
  const Config cfg(config);
  size_t ng = cfg.groups.ncol();
  if (ng == 0) {
    const Data d(data);
    if (d.Ny > 0) {
      ng = d.Ny;
    } else if (d.elgm_matrix.p.size() > 0) {
      ng = static_cast<size_t>(d.elgm_matrix.ncol());
    }
  }
  return ng > 0 ? ng : 1;
}

inline std::vector<GroupPack> build_ad_fun_obs(
  const Rcpp::List& data,
  const Rcpp::List& config,
  const std::string& obs_name) {

  LogDensObsFn log_dens = resolve_log_dens_obs(obs_name);
  const Config cfg(config);
  const size_t ng = count_obs_groups(data, config);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_obs " << obs_name << " groups " << ng << "\n";
  }

  std::vector<GroupPack> result(ng);

  CPPAD_TESTVECTOR(double) ad_params_G(cfg.Nparams);
  for (size_t d = 0; d < cfg.Nparams; ++d) {
    ad_params_G[d] = cfg.params[d];
  }

  CppAD::vector<CppAD::AD<double>> ad_params(cfg.Nparams);
  for (size_t d = 0; d < cfg.Nparams; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  for (size_t d = 0; d < ng; ++d) {
    CppAD::Independent(ad_params);
    auto result_here = log_dens(ad_params, data, config, d);
    CppAD::ADFun<double> fun(ad_params, result_here);
    result[d].fun = std::move(fun);
    adpack_sparsity(ad_params_G, cfg.Sgamma, result[d], cfg.verbose);
  }

  return result;
}

#endif
