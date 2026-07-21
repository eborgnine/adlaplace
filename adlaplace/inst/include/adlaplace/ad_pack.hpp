#ifndef ADLAPLACE_ADFUN_HPP
#define ADLAPLACE_ADFUN_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "adlaplace/ad_data.hpp"
#include "adlaplace/backend.hpp"
#include "adlaplace/rviews.hpp"

inline CPPAD_TESTVECTOR(double)
    make_ad_params_seed(const Config &cfg, const ad_data &model) {

  CPPAD_TESTVECTOR(double) ad_params_G(model.num_full);
  for (std::size_t d = 0; d < model.num_beta; ++d) {
    ad_params_G[d] = cfg.beta[d];
  }
  for (std::size_t d = 0; d < model.num_gamma; ++d) {
    ad_params_G[model.num_beta + d] = cfg.gamma[d];
  }
  for (std::size_t d = 0; d < model.num_theta; ++d) {
    ad_params_G[model.num_beta + model.num_gamma + d] = cfg.theta[d];
  }
  return ad_params_G;
}

static const std::string JAC_COLOR = "cppad";
static const std::string HESS_COLOR = "cppad.symmetric";

template <class SizeVec, class ValVec>
inline ValVec &rcv_val(CppAD::sparse_rcv<SizeVec, ValVec> &rcv) {
  return const_cast<ValVec &>(rcv.val());
}

inline const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &empty_sparse_rc() {
  static CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> empty;
  return empty;
}

inline void
adpack_discover_grad(GroupPack &gp,
                     CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &grad,
                     const std::size_t n_params) {

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> pattern_in;
  pattern_in.resize(n_params, n_params, n_params);
  for (std::size_t j = 0; j < n_params; ++j) {
    pattern_in.set(j, j, j);
  }

  const bool transpose = false;
  const bool dependency = false;
  const bool internal_bool = false;

  gp.fun.for_jac_sparsity(pattern_in, transpose, dependency, internal_bool,
                          grad);
}

inline void
adpack_discover_hessian(GroupPack &gp,
                        CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian,
                        const std::size_t n_params) {

  CppAD::vectorBool select_domain(n_params), select_range(1);
  for (std::size_t j = 0; j < n_params; ++j) {
    select_domain[j] = true;
  }
  select_range[0] = true;

  const bool internal_bool = false;
  gp.fun.for_hes_sparsity(select_domain, select_range, internal_bool, hessian);
}

inline void adpack_sparsity(const CPPAD_TESTVECTOR(double) & x,
                            const std::vector<int> &subset, GroupPack &gp,
                            const bool verbose = false,
                            const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
                                &hessian = empty_sparse_rc()) {

  const std::size_t n_params = x.size();

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad;
  adpack_discover_grad(gp, grad, n_params);

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_here;
  if (hessian.nnz() > 0) {
    if (verbose) {
      Rcpp::Rcout << "using external Hessian sparsity pattern\n";
    }
    hessian_here = hessian;
  } else {
    if (verbose) {
      Rcpp::Rcout << "discovering Hessian sparsity pattern\n";
    }
    adpack_discover_hessian(gp, hessian_here, n_params);
  }

  gp.w.resize(1);
  gp.w[0] = 1.0;
  gp.x.resize(n_params);

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad_inner;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_upper;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner_upper;

  const auto &full_cols = grad.col();

  std::vector<unsigned char> insubset(n_params, 0);
  for (std::size_t v : subset) {
    insubset[v] = 1;
  }

  {
    std::size_t k = 0;
    for (std::size_t j : full_cols) {
      k += insubset[j];
    }
    grad_inner.resize(1, n_params, k);

    std::size_t t = 0;
    for (std::size_t j : full_cols) {
      if (insubset[j])
        grad_inner.set(t++, 0, j);
    }
  }

  {
    const auto &r = hessian_here.row();
    const auto &c = hessian_here.col();
    const std::size_t k = hessian_here.nnz();

    std::size_t k_upper = 0, k_inner = 0, k_inner_upper = 0;

    for (std::size_t idx = 0; idx < k; ++idx) {
      const std::size_t i = r[idx], j = c[idx];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper)
        ++k_upper;
      if (is_inner)
        ++k_inner;
      if (is_upper && is_inner)
        ++k_inner_upper;
    }

    hessian_upper.resize(n_params, n_params, k_upper);
    hessian_inner.resize(n_params, n_params, k_inner);
    hessian_inner_upper.resize(n_params, n_params, k_inner_upper);

    std::size_t tu = 0, ti = 0, tiu = 0;
    for (std::size_t idx = 0; idx < k; ++idx) {
      const std::size_t i = r[idx], j = c[idx];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper)
        hessian_upper.set(tu++, i, j);
      if (is_inner)
        hessian_inner.set(ti++, i, j);
      if (is_upper && is_inner)
        hessian_inner_upper.set(tiu++, i, j);
    }
  }

  if (verbose) {
    Rcpp::Rcout << "  sparsity: grad " << grad.nnz() << ", grad_inner "
                << grad_inner.nnz() << ", hes full" << hessian_here.nnz()
                << ", hes upper" << hessian_upper.nnz() << ", hes_inner "
                << hessian_inner_upper.nnz() << "\n";
  }

  gp.fun.Forward(0, x);

  gp.pattern_grad =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          grad);
  gp.fun.sparse_jac_rev(x, gp.pattern_grad, grad, JAC_COLOR, gp.work_grad);

  gp.pattern_grad_inner =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          grad_inner);
  gp.fun.sparse_jac_rev(x, gp.pattern_grad_inner, grad_inner, JAC_COLOR,
                        gp.work_inner_grad);

  gp.pattern_hessian =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          hessian_upper);
  gp.fun.sparse_hes(x, gp.w, gp.pattern_hessian, hessian_here, HESS_COLOR,
                    gp.work_hess);

  gp.pattern_hessian_inner =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          hessian_inner_upper);
  gp.fun.sparse_hes(x, gp.w, gp.pattern_hessian_inner, hessian_inner,
                    HESS_COLOR, gp.work_inner_hess);
}

inline size_t count_obs_shards(const ad_data &model, const Rcpp::List &config) {
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

inline std::vector<GroupPack> build_ad_fun_obs(const ad_data &model,
                                               const Rcpp::List &config,
                                               LogDensObsFn log_dens) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model);
  const size_t ng = count_obs_shards(model, config);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_obs groups " << ng << "\n";
  }

  std::vector<GroupPack> result(ng);
  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  for (size_t d = 0; d < ng; ++d) {
    if (cfg.verbose) {
      Rcpp::Rcout << "  taping observation group " << (d + 1) << " / " << ng
                  << "\n";
    }
    CppAD::vector<CppAD::AD<double>> ad_params(model.num_full);
    for (size_t j = 0; j < model.num_full; ++j) {
      ad_params[j] = ad_params_G[j];
    }
    CppAD::Independent(ad_params);
    auto result_here = log_dens(ad_params, model, cfg, d);
    CppAD::ADFun<double> fun(ad_params, result_here);
    result[d].fun = std::move(fun);
  }

  for (size_t d = 0; d < ng; ++d) {
    if (cfg.verbose) {
      Rcpp::Rcout << "  sparsity observation group " << (d + 1) << " / " << ng
                  << "\n";
    }
    adpack_sparsity(ad_params_G, model.seq_gamma, result[d], cfg.verbose);
  }

  return result;
}

inline GroupPack build_ad_fun_parameters(
    const ad_data &model, const Rcpp::List &config, LogDensSingleDataFn log_dens,
    const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian = empty_sparse_rc()) {

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

  CppAD::vector<CppAD::AD<double>> result_here =
      log_dens(ad_params, model, cfg);

  CppAD::ADFun<double> fun(ad_params, result_here);

  GroupPack pack;
  pack.fun = std::move(fun);
  pack.owner_thread_assigned = false;
  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_parameters: computing sparsity...\n";
  }
  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose, hessian);
  return pack;
}

#endif
