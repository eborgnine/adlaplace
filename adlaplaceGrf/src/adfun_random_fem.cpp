#include <Rcpp.h>
#include <Rinternals.h>

#include "adlaplace/ad_data.hpp"
#include "adlaplace/adfun.hpp"
#include "adlaplace/adfun_random.hpp"
#include "adlaplace/atomics.hpp"
#include "adlaplace/chol_update_impl.hpp"
#include "adlaplace/extension.hpp"

#include <cmath>
#include <set>
#include <utility>
#include <vector>

// Declared by ADLAPLACE_DEFINE_BACKEND in backend.cpp
extern adlaplace_shard *adlaplace_grf_make_shard(GroupPack &&);

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

void validate_fem_alpha(const FemCholPayload &pay, const ad_data &model,
                        const int expect_alpha) {
  if (pay.alpha != expect_alpha) {
    Rcpp::stop("precision alpha (%d) does not match kernel alpha (%d)",
               pay.alpha, expect_alpha);
  }
  if (model.num_theta < 2) {
    Rcpp::stop("random_fem requires two theta parameters (range, sd)");
  }
}

struct FemThetaScale {
  std::size_t idx_range;
  std::size_t idx_sd;
  CppAD::AD<double> tau2;
  CppAD::AD<double> k2;
  CppAD::AD<double> k4;
  CppAD::AD<double> k6;
};

// Public theta is (range, sd) with the geostatsp/INLA practical range
//   range = rho = sqrt(8*nu)/kappa
// and sd = marginal field SD. Internally convert to SPDE (kappa, tau).
// In 2D, nu = alpha - 1 (alpha=2 => nu=1; alpha=3 => nu=2).
FemThetaScale fem_theta_scale(const CppAD::vector<CppAD::AD<double>> &x,
                              const ad_data &model, const Config &config,
                              const int expect_alpha) {
  FemThetaScale out;
  out.idx_range = model.theta_index(0);
  out.idx_sd = model.theta_index(1);
  CppAD::AD<double> log_range = x[out.idx_range];
  CppAD::AD<double> log_sd = x[out.idx_sd];
  if (!config.transform_theta) {
    log_range = CppAD::log(x[out.idx_range]);
    log_sd = CppAD::log(x[out.idx_sd]);
  }
  const CppAD::AD<double> range = CppAD::exp(log_range);
  const CppAD::AD<double> sd = CppAD::exp(log_sd);
  const double nu = static_cast<double>(expect_alpha - 1);
  const CppAD::AD<double> kappa =
      CppAD::sqrt(CppAD::AD<double>(8.0 * nu)) / range;
  out.k2 = kappa * kappa;
  out.k4 = out.k2 * out.k2;
  out.k6 = out.k4 * out.k2;

  // alpha=2 (nu=1): sigma^2 = 1/(4 pi kappa^2 tau^2)
  // alpha=3 (nu=2): sigma^2 = 1/(8 pi kappa^4 tau^2)
  CppAD::AD<double> tau;
  if (expect_alpha == 2) {
    tau = CppAD::AD<double>(1) /
          (kappa * sd * CppAD::sqrt(CppAD::AD<double>(4.0 * 3.14159265358979323846)));
  } else {
    tau = CppAD::AD<double>(1) /
          (out.k2 * sd * CppAD::sqrt(CppAD::AD<double>(8.0 * 3.14159265358979323846)));
  }
  out.tau2 = tau * tau;
  return out;
}

std::vector<CppAD::AD<double>>
assemble_Q_x(const FemCholPayload &pay, const FemThetaScale &th,
             const int expect_alpha) {
  const std::size_t nnz = pay.Q_i.size();
  std::vector<CppAD::AD<double>> Q_x(nnz);
  for (std::size_t k = 0; k < nnz; ++k) {
    if (expect_alpha == 2) {
      Q_x[k] = th.tau2 * (th.k4 * pay.C_x[k] +
                          CppAD::AD<double>(2) * th.k2 * pay.G_x[k] +
                          pay.G2_x[k]);
    } else {
      Q_x[k] = th.tau2 * (th.k6 * pay.C_x[k] +
                          CppAD::AD<double>(3) * th.k4 * pay.G_x[k] +
                          CppAD::AD<double>(3) * th.k2 * pay.G2_x[k] +
                          pay.G3_x[k]);
    }
  }
  return Q_x;
}

CppAD::vector<CppAD::AD<double>>
random_fem_ssq_impl(const CppAD::vector<CppAD::AD<double>> &x,
                    const ad_data &model, const Config &config,
                    const int expect_alpha) {
  const FemCholPayload pay = read_fem_payload(model);
  validate_fem_alpha(pay, model, expect_alpha);
  const FemThetaScale th = fem_theta_scale(x, model, config, expect_alpha);
  const std::vector<CppAD::AD<double>> Q_x = assemble_Q_x(pay, th, expect_alpha);

  const std::size_t n = pay.Q_p.size() - 1;
  const std::vector<std::size_t> gidx = model.all_gamma_global_indices();
  if (gidx.size() != n) {
    Rcpp::stop("random_fem: length(gamma) (%d) != nrow(Q) (%d)",
               static_cast<int>(gidx.size()), static_cast<int>(n));
  }

  CppAD::AD<double> qf = 0.0;
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
  result[0] = -CppAD::AD<double>(0.5) * qf;
  return result;
}

CppAD::vector<CppAD::AD<double>>
random_fem_det_impl(const CppAD::vector<CppAD::AD<double>> &x,
                    const ad_data &model, const Config &config,
                    const int expect_alpha) {
  const FemCholPayload pay = read_fem_payload(model);
  validate_fem_alpha(pay, model, expect_alpha);
  const FemThetaScale th = fem_theta_scale(x, model, config, expect_alpha);
  const std::vector<CppAD::AD<double>> Q_x = assemble_Q_x(pay, th, expect_alpha);

  const std::size_t n = pay.Q_p.size() - 1;
  std::vector<CppAD::AD<double>> L_x(pay.L1_i.size());
  std::vector<CppAD::AD<double>> D(n);
  const CppAD::AD<double> log_det = adlaplace::chol::chol_update_csc(
      pay.Q_p, pay.Q_i, Q_x, pay.perm, pay.L1_p, pay.L1_i, L_x, D);

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = CppAD::AD<double>(0.5) * log_det -
              CppAD::AD<double>(n * ONEHALFLOGTWOPI);
  return result;
}

CppAD::vector<CppAD::AD<double>>
random_fem_ssq_2(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
                 const Config &config) {
  return random_fem_ssq_impl(x, model, config, 2);
}

CppAD::vector<CppAD::AD<double>>
random_fem_ssq_3(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
                 const Config &config) {
  return random_fem_ssq_impl(x, model, config, 3);
}

CppAD::vector<CppAD::AD<double>>
random_fem_det_2(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
                 const Config &config) {
  return random_fem_det_impl(x, model, config, 2);
}

CppAD::vector<CppAD::AD<double>>
random_fem_det_3(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
                 const Config &config) {
  return random_fem_det_impl(x, model, config, 3);
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_ssq_sparsity(const ad_data &model) {
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

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_det_sparsity(const ad_data &model) {
  const std::size_t n_params = model.num_full;
  const std::size_t idx_kappa = model.theta_index(0);
  const std::size_t idx_tau = model.theta_index(1);

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  pairs.insert({idx_kappa, idx_kappa});
  pairs.insert({idx_tau, idx_tau});
  pairs.insert({idx_kappa, idx_tau});
  pairs.insert({idx_tau, idx_kappa});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, n_params, pairs);
  return hessian;
}

std::vector<GroupPack> build_random_fem_packs(const ad_data &ad_model,
                                              const Rcpp::List &config,
                                              LogDensSingleRandomFn ssq_fn,
                                              LogDensSingleRandomFn det_fn) {
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hes_ssq =
      random_fem_ssq_sparsity(ad_model);
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hes_det =
      random_fem_det_sparsity(ad_model);
  std::vector<GroupPack> packs;
  packs.reserve(2);
  packs.push_back(build_ad_fun_random(ad_model, config, ssq_fn, hes_ssq));
  packs.push_back(build_ad_fun_random(ad_model, config, det_fn, hes_det));
  return packs;
}

} // namespace

//' Build raw AD handle for a random_fem_2 term (ssq + det shards)
//'
//' @param model An \code{ad_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_fun_ptr} with two shards.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP create_ad_fun_random_fem_2(SEXP model, Rcpp::List config) {
  const ad_data ad_model(model);
  std::vector<GroupPack> packs =
      build_random_fem_packs(ad_model, config, random_fem_ssq_2, random_fem_det_2);
  ad_fun *groups = packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                   ad_model.num_theta, adlaplace_grf_make_shard);
  return make_ad_fun_ptr(groups);
}

//' Build raw AD handle for a random_fem_3 term (ssq + det shards)
//'
//' @param model An \code{ad_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_fun_ptr} with two shards.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP create_ad_fun_random_fem_3(SEXP model, Rcpp::List config) {
  const ad_data ad_model(model);
  std::vector<GroupPack> packs =
      build_random_fem_packs(ad_model, config, random_fem_ssq_3, random_fem_det_3);
  ad_fun *groups = packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                   ad_model.num_theta, adlaplace_grf_make_shard);
  return make_ad_fun_ptr(groups);
}
