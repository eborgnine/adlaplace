#include <Rcpp.h>
#include <Rinternals.h>

#include "adlaplace/density_data.hpp"
#include "adlaplace/ad_pack.hpp"
#include "adlaplace/ad_pack_random.hpp"
#include "adlaplace/chol_update_impl.hpp"
#include "adlaplace/extension.hpp"

#include <cmath>
#include <set>
#include <string>
#include <utility>
#include <vector>

// paste0("const double TWOPI = ", format(2 * Rmpfr::Const("pi", prec =
// 120)),";")
const double PIx4 = 12.566370614359172953850573533118011539;
const double PIx8 = 25.132741228718345907701147066236023078;
// paste0("const double ONEHALFLOGTWOPI = ", format(0.5 * log(2 *
//   Rmpfr::Const("pi", prec = 120))), ";")
const double ONEHALFLOGTWOPI = 0.91893853320467274178032973640561763976;

// Declared by ADLAPLACE_DEFINE_BACKEND in backend.cpp
extern ad_shard *adlaplace_fem_make_shard(AdTape &&);

namespace adlaplace {
namespace chol {
template <> struct CholLogDetTraits<CppAD::AD<double>> {
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

FemCholPayload read_fem_payload(const density_data &model) {
  if (Rf_isNull(model.precision) || TYPEOF(model.precision) != VECSXP) {
    Rcpp::stop(
        "random_fem precision must be a list from fem_precision_payload()");
  }
  Rcpp::List prec(model.precision);
  for (const char *key :
       {"Q_p", "Q_i", "C_x", "G_x", "G2_x", "chol", "alpha"}) {
    if (!prec.containsElementNamed(key)) {
      Rcpp::stop("random_fem precision missing '%s'", key);
    }
  }
  if (model.num_theta < 2) {
    Rcpp::stop("random_fem requires two theta parameters (range, sd)");
  }
  FemCholPayload out;
  out.alpha = Rcpp::as<int>(prec["alpha"]);
  out.Q_p = Rcpp::as<std::vector<int>>(prec["Q_p"]);
  out.Q_i = Rcpp::as<std::vector<int>>(prec["Q_i"]);
  out.C_x = Rcpp::as<std::vector<double>>(prec["C_x"]);
  out.G_x = Rcpp::as<std::vector<double>>(prec["G_x"]);
  out.G2_x = Rcpp::as<std::vector<double>>(prec["G2_x"]);
  out.G3_x = prec.containsElementNamed("G3_x")
                 ? Rcpp::as<std::vector<double>>(prec["G3_x"])
                 : std::vector<double>(out.Q_i.size(), 0.0);
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

// Public theta (range, sd): practical range rho = sqrt(8*nu)/kappa and field
// SD. In 2D, nu = Alpha - 1. Convert on tape to SPDE (kappa, tau).
template <int Alpha>
void fem_kappa_tau2(const CppAD::vector<CppAD::AD<double>> &x,
                    const density_data &model, const Config &config,
                    CppAD::AD<double> &k2, CppAD::AD<double> &k4,
                    CppAD::AD<double> &k6, CppAD::AD<double> &tau2) {
  const std::size_t t0_global = model.theta_index(0);
  const std::size_t t0_row = t0_global - model.num_beta - model.num_gamma;
  const std::size_t t1_global = model.theta_index(1);
  const std::size_t t1_row = t1_global - model.num_beta - model.num_gamma;
  CppAD::AD<double> log_range = x[t0_global];
  CppAD::AD<double> log_sd = x[t1_global];
  if (!transform_theta_at(config, t0_row)) {
    log_range = CppAD::log(log_range);
  }
  if (!transform_theta_at(config, t1_row)) {
    log_sd = CppAD::log(log_sd);
  }
  const CppAD::AD<double> range = CppAD::exp(log_range);
  const CppAD::AD<double> sd = CppAD::exp(log_sd);
  const CppAD::AD<double> kappa =
      CppAD::sqrt(CppAD::AD<double>(8.0 * (Alpha - 1))) / range;
  k2 = kappa * kappa;
  k4 = k2 * k2;
  k6 = k4 * k2;
  // Alpha=2: sigma^2 = 1/(4 pi kappa^2 tau^2)
  // Alpha=3: sigma^2 = 1/(8 pi kappa^4 tau^2)
  const CppAD::AD<double> tau =
      (Alpha == 2) ? CppAD::AD<double>(1) /
                         (kappa * sd * CppAD::sqrt(CppAD::AD<double>(PIx4)))
                   : CppAD::AD<double>(1) /
                         (k2 * sd * CppAD::sqrt(CppAD::AD<double>(PIx8)));
  tau2 = tau * tau;
}

template <int Alpha>
std::vector<CppAD::AD<double>>
assemble_Q_x(const FemCholPayload &pay, const CppAD::AD<double> &k2,
             const CppAD::AD<double> &k4, const CppAD::AD<double> &k6,
             const CppAD::AD<double> &tau2) {
  std::vector<CppAD::AD<double>> Q_x(pay.Q_i.size());
  for (std::size_t k = 0; k < Q_x.size(); ++k) {
    if constexpr (Alpha == 2) {
      Q_x[k] = tau2 * (k4 * pay.C_x[k] +
                       CppAD::AD<double>(2) * k2 * pay.G_x[k] + pay.G2_x[k]);
    } else {
      Q_x[k] =
          tau2 * (k6 * pay.C_x[k] + CppAD::AD<double>(3) * k4 * pay.G_x[k] +
                  CppAD::AD<double>(3) * k2 * pay.G2_x[k] + pay.G3_x[k]);
    }
  }
  return Q_x;
}

template <int Alpha>
CppAD::vector<CppAD::AD<double>>
random_fem_ssq(const CppAD::vector<CppAD::AD<double>> &x, const density_data &model,
               const Config &config) {
  const FemCholPayload pay = read_fem_payload(model);
  if (pay.alpha != Alpha) {
    Rcpp::stop("precision alpha (%d) does not match kernel alpha (%d)",
               pay.alpha, Alpha);
  }
  CppAD::AD<double> k2, k4, k6, tau2;
  fem_kappa_tau2<Alpha>(x, model, config, k2, k4, k6, tau2);
  const auto Q_x = assemble_Q_x<Alpha>(pay, k2, k4, k6, tau2);

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
      qf += (row == col ? CppAD::AD<double>(1) : CppAD::AD<double>(2)) *
            x[gidx[row]] * v * x[gidx[col]];
    }
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -CppAD::AD<double>(0.5) * qf;
  return result;
}

template <int Alpha>
CppAD::vector<CppAD::AD<double>>
random_fem_det(const CppAD::vector<CppAD::AD<double>> &x, const density_data &model,
               const Config &config) {
  const FemCholPayload pay = read_fem_payload(model);
  if (pay.alpha != Alpha) {
    Rcpp::stop("precision alpha (%d) does not match kernel alpha (%d)",
               pay.alpha, Alpha);
  }
  CppAD::AD<double> k2, k4, k6, tau2;
  fem_kappa_tau2<Alpha>(x, model, config, k2, k4, k6, tau2);
  const auto Q_x = assemble_Q_x<Alpha>(pay, k2, k4, k6, tau2);

  const std::size_t n = pay.Q_p.size() - 1;
  std::vector<CppAD::AD<double>> L_x(pay.L1_i.size());
  std::vector<CppAD::AD<double>> D(n);
  const CppAD::AD<double> log_det = adlaplace::chol::chol_update_csc(
      pay.Q_p, pay.Q_i, Q_x, pay.perm, pay.L1_p, pay.L1_i, L_x, D);

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] =
      CppAD::AD<double>(0.5) * log_det - CppAD::AD<double>(n * ONEHALFLOGTWOPI);
  return result;
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_ssq_sparsity(const density_data &model) {
  const FemCholPayload pay = read_fem_payload(model);
  const std::size_t idx_range = model.theta_index(0);
  const std::size_t idx_sd = model.theta_index(1);
  const std::vector<std::size_t> gidx = model.all_gamma_global_indices();

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  for (std::size_t col = 0; col < gidx.size(); ++col) {
    for (int pos = pay.Q_p[col]; pos < pay.Q_p[col + 1]; ++pos) {
      const std::size_t row =
          static_cast<std::size_t>(pay.Q_i[static_cast<std::size_t>(pos)]);
      pairs.insert({gidx[row], gidx[col]});
      pairs.insert({gidx[col], gidx[row]});
    }
  }
  for (std::size_t k = 0; k < gidx.size(); ++k) {
    pairs.insert({gidx[k], idx_range});
    pairs.insert({idx_range, gidx[k]});
    pairs.insert({gidx[k], idx_sd});
    pairs.insert({idx_sd, gidx[k]});
  }
  pairs.insert({idx_range, idx_range});
  pairs.insert({idx_sd, idx_sd});
  pairs.insert({idx_range, idx_sd});
  pairs.insert({idx_sd, idx_range});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, model.num_full, pairs);
  return hessian;
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_det_sparsity(const density_data &model) {
  const std::size_t idx_range = model.theta_index(0);
  const std::size_t idx_sd = model.theta_index(1);

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  pairs.insert({idx_range, idx_range});
  pairs.insert({idx_sd, idx_sd});
  pairs.insert({idx_range, idx_sd});
  pairs.insert({idx_sd, idx_range});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, model.num_full, pairs);
  return hessian;
}

template <int Alpha>
SEXP create_ad_shard_random_fem_ssq(SEXP model, Rcpp::List config) {
  const density_data ad_model(model);
  std::vector<AdTape> packs;
  packs.push_back(build_ad_fun_random(ad_model, config, random_fem_ssq<Alpha>,
                                      random_fem_ssq_sparsity(ad_model)));
  return make_ad_pack_ptr(packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                         ad_model.num_theta,
                                         adlaplace_fem_make_shard));
}

LogDensSingleDataFn resolve_fem_det(const std::string &name) {
  if (name == "random_fem_det_2") {
    return random_fem_det<2>;
  }
  if (name == "random_fem_det_3") {
    return random_fem_det<3>;
  }
  Rcpp::stop("unknown parameters density: %s", name.c_str());
}

} // namespace

//' Build raw AD handle for a random_fem_ssq_2 term
//'
//' @param model An \code{density_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP create_ad_shard_random_fem_ssq_2(SEXP model, Rcpp::List config) {
  return create_ad_shard_random_fem_ssq<2>(model, config);
}

//' Build raw AD handle for a random_fem_ssq_3 term
//'
//' @param model An \code{density_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP create_ad_shard_random_fem_ssq_3(SEXP model, Rcpp::List config) {
  return create_ad_shard_random_fem_ssq<3>(model, config);
}

//' Build parameters-shard \code{ad_pack_ptr} for FEM log-determinant densities.
//'
//' @param model \code{density_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @param name Density name (\code{"random_fem_det_2"} or \code{"random_fem_det_3"}).
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP get_ad_pack_raw_parameters(SEXP model, Rcpp::List config, std::string name) {
  const density_data ad_model(model);
  AdTape pack = build_ad_fun_parameters(
      ad_model, config, resolve_fem_det(name),
      random_fem_det_sparsity(ad_model));
  std::vector<AdTape> packs;
  packs.push_back(std::move(pack));
  return make_ad_pack_ptr(packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                         ad_model.num_theta,
                                         adlaplace_fem_make_shard));
}
