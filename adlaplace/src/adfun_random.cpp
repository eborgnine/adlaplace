#include <Rcpp.h>
#include <Rinternals.h>

#include "adlaplace/ad_data.hpp"
#include "adlaplace/adfun.hpp"
#include "adlaplace/atomics.hpp"
#include "adlaplace/register.hpp"
#include "adlaplace/rviews.hpp"

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

extern adlaplace_shard *adlaplace_make_shard(GroupPack &&);

namespace {

CppAD::vector<CppAD::AD<double>>
random_diagonal_impl(const CppAD::vector<CppAD::AD<double>> &x,
                     const ad_data &model, const NumVecView &Q,
                     const std::vector<std::size_t> &gamma_indices,
                     const Config &config) {

  const std::size_t Ngamma = gamma_indices.size();
  if (Q.size() != Ngamma) {
    Rcpp::warning("precision length (%d) differs from gamma_map columns (%d)",
                  static_cast<int>(Q.size()), static_cast<int>(Ngamma));
  }

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> logSd;
  if (config.transform_theta) {
    logSd = x[theta_index];
  } else {
    logSd = CppAD::log(x[theta_index]);
  }

  CppAD::AD<double> precision = CppAD::exp(-2 * logSd);
  CppAD::AD<double> qpart = 0.0;
  double logQsum = 0.0;

  const std::size_t nq = std::min(Ngamma, Q.size());
  for (std::size_t k = 0; k < nq; ++k) {
    const std::size_t gidx = gamma_indices[k];
    const double qk = Q[k];
    qpart += x[gidx] * x[gidx] * qk;
    logQsum += std::log(qk);
  }
  qpart *= CppAD::AD<double>(0.5) * precision;

  CppAD::AD<double> qDet = logSd * CppAD::AD<double>(Ngamma) +
                           CppAD::AD<double>(Ngamma * ONEHALFLOGTWOPI);

  if (config.verbose) {
    Rcpp::Rcout << "theta index " << theta_index << " logVariance " << logSd
                << " precision " << precision << "\n";
    Rcpp::Rcout << "random_diagonal n_gamma " << Ngamma << " qDet " << qDet
                << " qpart " << qpart << "\n";
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -qpart - qDet + CppAD::AD<double>(0.5 * logQsum);
  return result;
}

// CppAD's sparse_hes hands the pattern directly to the cppad.symmetric
// coloring, which treats missing entries as structural zeros; it does NOT
// symmetrize. Analytic patterns must therefore contain BOTH triangles,
// like the pattern for_hes_sparsity discovers. adpack_sparsity extracts
// the upper-triangle subsets itself.
void set_sparse_rc_pairs(
    CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian,
    const std::size_t n_params,
    const std::set<std::pair<std::size_t, std::size_t>> &pairs) {
  hessian.resize(n_params, n_params, pairs.size());
  std::size_t k = 0;
  for (const auto &rc : pairs) {
    hessian.set(k++, rc.first, rc.second);
  }
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_diagonal_sparsity(const ad_data &model) {

  const std::vector<std::size_t> gamma_indices =
      model.all_gamma_global_indices();
  const std::size_t theta_index = model.theta_index(0);
  const std::size_t n_gamma = gamma_indices.size();
  const std::size_t n_params = model.num_full;

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  hessian.resize(n_params, n_params, 3 * n_gamma + 1);
  hessian.set(0, theta_index, theta_index);
  std::size_t t = 1;
  for (std::size_t k = 0; k < n_gamma; ++k) {
    const std::size_t gi = gamma_indices[k];
    hessian.set(t++, gi, gi);
    hessian.set(t++, gi, theta_index);
    hessian.set(t++, theta_index, gi);
  }
  return hessian;
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_mult_sparsity(const ad_data &model) {

  const DgCView Q = model.mult_precision_Q();
  const std::size_t theta_index = model.theta_index(0);
  const std::size_t n_params = model.num_full;
  const int n_term = model.gamma_map.ncol();

  std::set<std::pair<std::size_t, std::size_t>> pairs;

  // gamma-gamma block: structural nonzeros of Q, mapped to global indices.
  // Q uses general (full symmetric) storage, but insert both orders anyway
  // in case only one triangle is stored.
  for (int col = 0; col < n_term; ++col) {
    const std::vector<std::size_t> gj_idx = model.gamma_global_indices(col);
    if (gj_idx.empty()) {
      continue;
    }
    const std::size_t gj = gj_idx[0];

    for (int k = Q.p[col]; k < Q.p[col + 1]; ++k) {
      const std::vector<std::size_t> gi_idx =
          model.gamma_global_indices(static_cast<int>(Q.i[k]));
      if (gi_idx.empty()) {
        continue;
      }
      const std::size_t gi = gi_idx[0];
      pairs.insert({gi, gj});
      pairs.insert({gj, gi});
    }
  }

  // gamma-theta cross terms (both orders) and the theta diagonal.
  for (int col = 0; col < n_term; ++col) {
    const std::vector<std::size_t> g_idx = model.gamma_global_indices(col);
    if (g_idx.empty()) {
      continue;
    }
    pairs.insert({g_idx[0], theta_index});
    pairs.insert({theta_index, g_idx[0]});
  }
  pairs.insert({theta_index, theta_index});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, n_params, pairs);
  return hessian;
}

GroupPack build_ad_fun_random(const ad_data &model, const Rcpp::List &config,
                              LogDensSingleRandomFn log_dens,
                              const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
                                  &hessian = empty_sparse_rc()) {

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

  GroupPack pack;
  pack.fun = std::move(fun);
  pack.owner_thread_assigned = false;
  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_random: computing sparsity...\n";
  }

  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose, hessian);
  return pack;
}

} // namespace

CppAD::vector<CppAD::AD<double>>
random_diagonal(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
                const Config &config) {

  if (Rf_isNull(model.precision)) {
    Rcpp::stop("precision is required for random_diagonal");
  }
  const NumVecView Q(model.precision);
  const std::vector<std::size_t> gamma_indices =
      model.all_gamma_global_indices();
  return random_diagonal_impl(x, model, Q, gamma_indices, config);
}

CppAD::vector<CppAD::AD<double>>
random_mult(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
            const Config &config) {

  const DgCView Q = model.mult_precision_Q();
  const int n_term = model.gamma_map.ncol();
  if (static_cast<std::size_t>(Q.nrow()) != static_cast<std::size_t>(n_term) ||
      static_cast<std::size_t>(Q.ncol()) != static_cast<std::size_t>(n_term)) {
    Rcpp::stop("random_mult Q is %d x %d but gamma_map has %d columns",
               Q.nrow(), Q.ncol(), n_term);
  }

  const double rank = model.mult_precision_rank();
  const double log_det = model.mult_precision_log_det();

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> logSd;
  if (config.transform_theta) {
    logSd = x[theta_index];
  } else {
    logSd = CppAD::log(x[theta_index]);
  }
  const CppAD::AD<double> tau = CppAD::exp(-2 * logSd);

  CppAD::AD<double> qf = 0.0;
  for (int j = 0; j < Q.ncol(); ++j) {
    const std::vector<std::size_t> gj_idx = model.gamma_global_indices(j);
    if (gj_idx.empty()) {
      Rcpp::stop("gamma_map column %d has no structural nonzero", j + 1);
    }
    const std::size_t gj = gj_idx[0];
    const int p0 = Q.p[j];
    const int p1 = Q.p[j + 1];
    for (int k = p0; k < p1; ++k) {
      const std::vector<std::size_t> gi_idx =
          model.gamma_global_indices(static_cast<int>(Q.i[k]));
      if (gi_idx.empty()) {
        Rcpp::stop("gamma_map column %d has no structural nonzero", Q.i[k] + 1);
      }
      const std::size_t gi = gi_idx[0];
      qf += x[gi] * Q.value(k) * x[gj];
    }
  }
  const CppAD::AD<double> qpart = CppAD::AD<double>(0.5) * tau * qf;

  if (config.verbose) {
    Rcpp::Rcout << "random_mult n_gamma " << n_term << " rank " << rank
                << " log_det " << log_det << " theta index " << theta_index
                << "\n";
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -qpart - logSd * CppAD::AD<double>(rank) +
              CppAD::AD<double>(0.5 * log_det - rank * ONEHALFLOGTWOPI);
  return result;
}

//' Build raw AD handle for a random_diagonal shard
//'
//' @param model An \code{ad_data} S4 object with \code{precision} slot set.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP create_ad_fun_random_diagonal(SEXP model, Rcpp::List config) {
  const ad_data ad_model(model);
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hes_pat =
      random_diagonal_sparsity(ad_model);
  GroupPack pack =
      build_ad_fun_random(ad_model, config, random_diagonal, hes_pat);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  ad_fun *groups = packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                   ad_model.num_theta, adlaplace_make_shard);
  return make_ad_fun_ptr(groups);
}

//' Build raw AD handle for a random_mult shard
//'
//' @param model An \code{ad_data} S4 object with \code{precision} slot set.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP create_ad_fun_random_mult(SEXP model, Rcpp::List config) {
  const ad_data ad_model(model);
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hes_pat =
      random_mult_sparsity(ad_model);
  GroupPack pack =
      build_ad_fun_random(ad_model, config, random_mult, hes_pat);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  ad_fun *groups = packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                   ad_model.num_theta, adlaplace_make_shard);
  return make_ad_fun_ptr(groups);
}

#include "adlaplace/chol_update_impl.hpp"

namespace adlaplace {
namespace chol {
template <>
struct CholLogDetTraits<CppAD::AD<double>> {
  static bool invalid_diag(const CppAD::AD<double> &) { return false; }
  static CppAD::AD<double> invalid_result() { return CppAD::AD<double>(0); }
  static CppAD::AD<double> log_of(const CppAD::AD<double> &d) {
    return CppAD::log(d);
  }
};
} // namespace chol
} // namespace adlaplace

namespace {

struct FemCholPayload {
  std::vector<int> Q_p;
  std::vector<int> Q_i;
  std::vector<double> C_x;
  std::vector<double> G_x;
  std::vector<double> G2_x;
  std::vector<double> G3_x;
  std::vector<int> perm;
  std::vector<int> L1_p;
  std::vector<int> L1_i;
  int alpha = 2;
};

FemCholPayload read_fem_payload(const ad_data &model) {
  if (Rf_isNull(model.precision) || TYPEOF(model.precision) != VECSXP) {
    Rcpp::stop("random_fem precision must be a list from fem_precision_payload()");
  }
  Rcpp::List prec(model.precision);
  for (const char *key : {"Q_p", "Q_i", "C_x", "G_x", "G2_x", "chol", "alpha"}) {
    if (!prec.containsElementNamed(key)) {
      Rcpp::stop("random_fem precision missing '%s'", key);
    }
  }
  FemCholPayload out;
  out.alpha = Rcpp::as<int>(prec["alpha"]);
  out.Q_p = Rcpp::as<std::vector<int>>(prec["Q_p"]);
  out.Q_i = Rcpp::as<std::vector<int>>(prec["Q_i"]);
  out.C_x = Rcpp::as<std::vector<double>>(prec["C_x"]);
  out.G_x = Rcpp::as<std::vector<double>>(prec["G_x"]);
  out.G2_x = Rcpp::as<std::vector<double>>(prec["G2_x"]);
  if (prec.containsElementNamed("G3_x")) {
    out.G3_x = Rcpp::as<std::vector<double>>(prec["G3_x"]);
  } else {
    out.G3_x.assign(out.Q_i.size(), 0.0);
  }
  Rcpp::List chol = prec["chol"];
  if (!chol.containsElementNamed("perm") || !chol.containsElementNamed("L1")) {
    Rcpp::stop("random_fem chol must contain perm and L1");
  }
  out.perm = Rcpp::as<std::vector<int>>(chol["perm"]);
  Rcpp::S4 L1 = chol["L1"];
  Rcpp::IntegerVector Lp = L1.slot("p");
  Rcpp::IntegerVector Li = L1.slot("i");
  out.L1_p.assign(Lp.begin(), Lp.end());
  out.L1_i.assign(Li.begin(), Li.end());
  if (out.C_x.size() != out.Q_i.size() || out.G_x.size() != out.Q_i.size() ||
      out.G2_x.size() != out.Q_i.size() || out.G3_x.size() != out.Q_i.size()) {
    Rcpp::stop("random_fem Gram coefficient vectors must match length(Q_i)");
  }
  return out;
}

CppAD::vector<CppAD::AD<double>>
random_fem_impl(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
                const Config &config, const int expect_alpha) {
  const FemCholPayload pay = read_fem_payload(model);
  if (pay.alpha != expect_alpha) {
    Rcpp::stop("precision alpha (%d) does not match kernel alpha (%d)",
               pay.alpha, expect_alpha);
  }
  if (model.num_theta < 2) {
    Rcpp::stop("random_fem requires two theta parameters (log_kappa, log_tau)");
  }

  const std::size_t idx_kappa = model.theta_index(0);
  const std::size_t idx_tau = model.theta_index(1);
  CppAD::AD<double> log_kappa = x[idx_kappa];
  CppAD::AD<double> log_tau = x[idx_tau];
  if (!config.transform_theta) {
    log_kappa = CppAD::log(x[idx_kappa]);
    log_tau = CppAD::log(x[idx_tau]);
  }
  const CppAD::AD<double> kappa = CppAD::exp(log_kappa);
  const CppAD::AD<double> tau = CppAD::exp(log_tau);
  const CppAD::AD<double> tau2 = tau * tau;
  const CppAD::AD<double> k2 = kappa * kappa;
  const CppAD::AD<double> k4 = k2 * k2;
  const CppAD::AD<double> k6 = k4 * k2;

  const std::size_t nnz = pay.Q_i.size();
  std::vector<CppAD::AD<double>> Q_x(nnz);
  for (std::size_t k = 0; k < nnz; ++k) {
    if (expect_alpha == 2) {
      Q_x[k] = tau2 * (k4 * pay.C_x[k] + CppAD::AD<double>(2) * k2 * pay.G_x[k] +
                       pay.G2_x[k]);
    } else {
      Q_x[k] = tau2 * (k6 * pay.C_x[k] + CppAD::AD<double>(3) * k4 * pay.G_x[k] +
                       CppAD::AD<double>(3) * k2 * pay.G2_x[k] + pay.G3_x[k]);
    }
  }

  const std::size_t n = pay.Q_p.size() - 1;
  std::vector<CppAD::AD<double>> L_x(pay.L1_i.size());
  std::vector<CppAD::AD<double>> D(n);
  const CppAD::AD<double> log_det = adlaplace::chol::chol_update_csc(
      pay.Q_p, pay.Q_i, Q_x, pay.perm, pay.L1_p, pay.L1_i, L_x, D);

  CppAD::AD<double> qf = 0.0;
  const std::vector<std::size_t> gidx = model.all_gamma_global_indices();
  if (gidx.size() != n) {
    Rcpp::stop("random_fem: length(gamma) (%d) != nrow(Q) (%d)",
               static_cast<int>(gidx.size()), static_cast<int>(n));
  }
  for (std::size_t col = 0; col < n; ++col) {
    for (int pos = pay.Q_p[col]; pos < pay.Q_p[col + 1]; ++pos) {
      const std::size_t row =
          static_cast<std::size_t>(pay.Q_i[static_cast<std::size_t>(pos)]);
      const CppAD::AD<double> v = Q_x[static_cast<std::size_t>(pos)];
      if (row == col) {
        qf += x[gidx[row]] * v * x[gidx[col]];
      } else {
        qf += CppAD::AD<double>(2) * x[gidx[row]] * v * x[gidx[col]];
      }
    }
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = CppAD::AD<double>(0.5) * log_det - CppAD::AD<double>(0.5) * qf -
              CppAD::AD<double>(n * ONEHALFLOGTWOPI);
  return result;
}

CppAD::vector<CppAD::AD<double>>
random_fem_2(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
             const Config &config) {
  return random_fem_impl(x, model, config, 2);
}

CppAD::vector<CppAD::AD<double>>
random_fem_3(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
             const Config &config) {
  return random_fem_impl(x, model, config, 3);
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_sparsity(const ad_data &model) {
  const FemCholPayload pay = read_fem_payload(model);
  const std::size_t n_params = model.num_full;
  const std::size_t idx_kappa = model.theta_index(0);
  const std::size_t idx_tau = model.theta_index(1);
  const std::vector<std::size_t> gidx = model.all_gamma_global_indices();
  const std::size_t n = gidx.size();

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  for (std::size_t col = 0; col < n; ++col) {
    for (int pos = pay.Q_p[col]; pos < pay.Q_p[col + 1]; ++pos) {
      const std::size_t row =
          static_cast<std::size_t>(pay.Q_i[static_cast<std::size_t>(pos)]);
      pairs.insert({gidx[row], gidx[col]});
      pairs.insert({gidx[col], gidx[row]});
    }
  }
  for (std::size_t k = 0; k < n; ++k) {
    pairs.insert({gidx[k], idx_kappa});
    pairs.insert({idx_kappa, gidx[k]});
    pairs.insert({gidx[k], idx_tau});
    pairs.insert({idx_tau, gidx[k]});
  }
  pairs.insert({idx_kappa, idx_kappa});
  pairs.insert({idx_tau, idx_tau});
  pairs.insert({idx_kappa, idx_tau});
  pairs.insert({idx_tau, idx_kappa});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, n_params, pairs);
  return hessian;
}

} // namespace

//' Build raw AD handle for a random_fem_2 shard
//'
//' @param model An \code{ad_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP create_ad_fun_random_fem_2(SEXP model, Rcpp::List config) {
  const ad_data ad_model(model);
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hes_pat =
      random_fem_sparsity(ad_model);
  GroupPack pack =
      build_ad_fun_random(ad_model, config, random_fem_2, hes_pat);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  ad_fun *groups = packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                   ad_model.num_theta, adlaplace_make_shard);
  return make_ad_fun_ptr(groups);
}

//' Build raw AD handle for a random_fem_3 shard
//'
//' @param model An \code{ad_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP create_ad_fun_random_fem_3(SEXP model, Rcpp::List config) {
  const ad_data ad_model(model);
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hes_pat =
      random_fem_sparsity(ad_model);
  GroupPack pack =
      build_ad_fun_random(ad_model, config, random_fem_3, hes_pat);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  ad_fun *groups = packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                   ad_model.num_theta, adlaplace_make_shard);
  return make_ad_fun_ptr(groups);
}
