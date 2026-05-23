#ifndef ADFUN_SINGLE_HPP
#define ADFUN_SINGLE_HPP

#include "adlaplace/api/density_registry.hpp"
#include "adlaplace/creators/adfun_common.hpp"
#include "adlaplace/creators/ad_model.hpp"

inline GroupPack build_ad_fun_random(
  const ad_model& model,
  const Rcpp::List& precision,
  const Rcpp::List& config,
  const std::string& single_name) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model, false);
  if (precision.containsElementNamed("Q")) {
    const NumVecView Q(precision["Q"]);
    const std::vector<std::size_t> gidx = model.all_gamma_global_indices();
    if (Q.size() != static_cast<R_xlen_t>(gidx.size())) {
      Rcpp::stop(
        "length(precision$Q) (%d) must match active gamma_map entries (%d)",
        static_cast<int>(Q.size()),
        static_cast<int>(gidx.size())
      );
    }
  }

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_random " << single_name << "\n";
  }

  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  CppAD::vector<CppAD::AD<double>> ad_params(model.num_full);
  for (size_t d = 0; d < model.num_full; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  CppAD::Independent(ad_params);

  LogDensSingleRandomDiagFn log_dens = resolve_log_dens_single_random_diag(single_name);
  CppAD::vector<CppAD::AD<double>> result_here = log_dens(ad_params, model, precision, config);

  CppAD::ADFun<double> fun(ad_params, result_here);

  GroupPack pack;
  pack.fun = std::move(fun);
  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose);
  return pack;
}

inline GroupPack build_ad_fun_parameters(
  const ad_model& model,
  const Rcpp::List& config,
  const std::string& single_name) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_parameters " << single_name << "\n";
  }

  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  CppAD::vector<CppAD::AD<double>> ad_params(model.num_full);
  for (size_t d = 0; d < model.num_full; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  CppAD::Independent(ad_params);

  LogDensSingleDataFn log_dens = resolve_log_dens_single_data(single_name);
  CppAD::vector<CppAD::AD<double>> result_here = log_dens(ad_params, model, config);

  CppAD::ADFun<double> fun(ad_params, result_here);

  GroupPack pack;
  pack.fun = std::move(fun);
  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose);
  return pack;
}

inline GroupPack build_ad_fun_single(
  const ad_model& model,
  const Rcpp::List& config,
  const std::string& single_name) {

  if (log_dens_single_uses_random_diag(single_name)) {
    Rcpp::stop("build_ad_fun_single requires precision for random density '%s'", single_name.c_str());
  }
  return build_ad_fun_parameters(model, config, single_name);
}

#endif
