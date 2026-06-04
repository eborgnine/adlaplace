#ifndef ADFUN_OBS_HPP
#define ADFUN_OBS_HPP

#include "adlaplace/api/density_registry.hpp"
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
  const std::string& obs_name) {

  LogDensObsFn log_dens = resolve_log_dens_obs(obs_name);
  const Config cfg(config);
  validate_config_matches_model(cfg, model);
  const size_t ng = count_obs_shards(model, config);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_obs " << obs_name << " groups " << ng << "\n";
  }

  std::vector<GroupPack> result(ng);
  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

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

  const int num_threads = cfg.num_threads > 0 ? cfg.num_threads : 1;
  if (num_threads > 1 && ng > 1) {
    cppad_parallel_setup(static_cast<std::size_t>(num_threads));
#pragma omp parallel for num_threads(num_threads) schedule(static,1)
    for (int di = 0; di < static_cast<int>(ng); ++di) {
      const size_t d = static_cast<size_t>(di);
      result[d].owner_thread = static_cast<std::size_t>(di % num_threads);
      adpack_sparsity(ad_params_G, model.seq_gamma, result[d], cfg.verbose);
    }
  } else {
    for (size_t d = 0; d < ng; ++d) {
      result[d].owner_thread = 0;
      adpack_sparsity(ad_params_G, model.seq_gamma, result[d], cfg.verbose);
    }
  }

  return result;
}

#endif
