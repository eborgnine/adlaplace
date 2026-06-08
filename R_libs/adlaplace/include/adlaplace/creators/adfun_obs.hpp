#ifndef ADFUN_OBS_HPP
#define ADFUN_OBS_HPP

#include "adlaplace/api/log_dens_fn.hpp"
#include "adlaplace/creators/adfun_common.hpp"
#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/ompad.hpp"

inline size_t count_obs_shards(const ad_data& model, const Rcpp::List& config) {
  const Config cfg(config);
  size_t ng = cfg.shards.ncol();
  if (ng == 0) {
    if (model.y.size() > 0) {
      ng = 1;
    } else if (model.elgm_matrix.p.size() > 0) {
      ng = static_cast<size_t>(model.elgm_matrix.ncol());
    }
  }
  return ng > 0 ? ng : 1;
}

inline std::vector<GroupPack> build_ad_fun_obs(
  const ad_data& model,
  const Rcpp::List& config,
  LogDensObsFn log_dens) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model);
  const size_t ng = count_obs_shards(model, config);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_obs groups " << ng << "\n";
  }

  std::vector<GroupPack> result(ng);
  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  cppad_parallel_setup(1);

  for (size_t d = 0; d < ng; ++d) {
    CppAD::vector<CppAD::AD<double>> ad_params(model.num_full);
    for (size_t j = 0; j < model.num_full; ++j) {
      ad_params[j] = ad_params_G[j];
    }
    CppAD::Independent(ad_params);
    auto result_here = log_dens(ad_params, model, config, d);
    CppAD::ADFun<double> fun(ad_params, result_here);
    result[d].fun = std::move(fun);
  }

  for (size_t d = 0; d < ng; ++d) {
    adpack_sparsity(ad_params_G, model.seq_gamma, result[d], cfg.verbose);
  }
  cppad_parallel_setup(1);

  return result;
}

#endif
