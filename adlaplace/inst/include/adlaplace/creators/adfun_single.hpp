#ifndef ADFUN_SINGLE_HPP
#define ADFUN_SINGLE_HPP

#include "adlaplace/api/log_dens_fn.hpp"
#include "adlaplace/creators/adfun_common.hpp"
#include "adlaplace/creators/ad_data.hpp"

inline GroupPack build_ad_fun_random(
  const ad_data& model,
  SEXP precision,
  const Rcpp::List& config,
  LogDensSingleRandomDiagFn log_dens) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model, false);
  if (Rf_isNull(precision)) {
    Rcpp::stop("precision is required for random densities");
  }
  if (TYPEOF(precision) == REALSXP || TYPEOF(precision) == INTSXP) {
    const NumVecView Q(precision);
    const int n_gamma_cols = model.gamma_map.ncol();
    if (Q.size() != static_cast<R_xlen_t>(n_gamma_cols)) {
      Rcpp::stop(
        "length(precision) (%d) must match ncol(gamma_map) (%d)",
        static_cast<int>(Q.size()),
        n_gamma_cols
      );
    }
  }

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_random: taping...\n";
  }

  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  CppAD::vector<CppAD::AD<double>> ad_params(model.num_full);
  for (size_t d = 0; d < model.num_full; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  CppAD::Independent(ad_params);

  CppAD::vector<CppAD::AD<double>> result_here = log_dens(ad_params, model, precision, config);

  CppAD::ADFun<double> fun(ad_params, result_here);

  GroupPack pack;
  pack.fun = std::move(fun);
  pack.owner_thread_assigned = false;
  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_random: computing sparsity...\n";
  }
  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose);
  return pack;
}

inline GroupPack build_ad_fun_parameters(
  const ad_data& model,
  const Rcpp::List& config,
  LogDensSingleDataFn log_dens) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_parameters: taping...\n";
  }

  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  CppAD::vector<CppAD::AD<double>> ad_params(model.num_full);
  for (size_t d = 0; d < model.num_full; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  CppAD::Independent(ad_params);

  CppAD::vector<CppAD::AD<double>> result_here = log_dens(ad_params, model, config);

  CppAD::ADFun<double> fun(ad_params, result_here);

  GroupPack pack;
  pack.fun = std::move(fun);
  pack.owner_thread_assigned = false;
  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_parameters: computing sparsity...\n";
  }
  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose);
  return pack;
}

#endif
