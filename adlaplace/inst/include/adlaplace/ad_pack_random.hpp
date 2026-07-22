#ifndef ADLAPLACE_ADFUN_RANDOM_HPP
#define ADLAPLACE_ADFUN_RANDOM_HPP

// Public helpers for random_* AD fun builders (diagonal / mult / FEM backends).

#include "adlaplace/ad_pack.hpp"
#include "adlaplace/backend.hpp"
#include "adlaplace/rviews.hpp"

#include <Rinternals.h>

#include <cstddef>
#include <set>
#include <utility>

// CppAD's sparse_hes hands the pattern directly to the cppad.symmetric
// coloring, which treats missing entries as structural zeros; it does NOT
// symmetrize. Analytic patterns must therefore contain BOTH triangles,
// like the pattern for_hes_sparsity discovers. adpack_sparsity extracts
// the upper-triangle subsets itself.
inline void set_sparse_rc_pairs(
    CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian,
    const std::size_t n_params,
    const std::set<std::pair<std::size_t, std::size_t>> &pairs) {
  hessian.resize(n_params, n_params, pairs.size());
  std::size_t k = 0;
  for (const auto &rc : pairs) {
    hessian.set(k++, rc.first, rc.second);
  }
}

inline AdTape build_ad_fun_random(
    const density_data &model, const Rcpp::List &config, LogDensSingleRandomFn log_dens,
    const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian = empty_sparse_rc()) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model, false);
  if (Rf_isNull(model.precision)) {
    Rcpp::stop("precision is required for random densities");
  }
  if (TYPEOF(model.precision) == REALSXP || TYPEOF(model.precision) == INTSXP) {
    const NumVecView Q(model.precision);
    const int n_gamma_cols = model.gamma_map.ncol();
    if (Q.size() != static_cast<R_xlen_t>(n_gamma_cols)) {
      Rcpp::stop("length(precision) (%d) must match ncol(gamma_map) (%d)",
                 static_cast<int>(Q.size()), n_gamma_cols);
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

  CppAD::vector<CppAD::AD<double>> result_here =
      log_dens(ad_params, model, cfg);

  CppAD::ADFun<double> fun(ad_params, result_here);

  AdTape pack;
  pack.fun = std::move(fun);
  pack.owner_thread_assigned = false;
  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_random: computing sparsity...\n";
  }

  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose, hessian);
  return pack;
}

#endif
