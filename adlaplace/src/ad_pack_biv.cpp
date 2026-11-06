#include <Rcpp.h>
#include <Rinternals.h>

#include "adlaplace/ad_pack.hpp"
#include "adlaplace/ad_pack_random.hpp"
#include "adlaplace/atomics.hpp"
#include "adlaplace/density_data.hpp"
#include "adlaplace/register.hpp"

#include <cmath>
#include <set>
#include <utility>

extern ad_shard* adlaplace_make_ad_shard(AdTape&&);

namespace {

int biv_n_pollutant(const density_data& model) {
  int n_pol = 2;
  if (TYPEOF(model.precision) == VECSXP) {
    Rcpp::List prec(model.precision);
    if (prec.containsElementNamed("n_pollutant")) {
      n_pol = Rcpp::as<int>(prec["n_pollutant"]);
    }
  }
  if (n_pol != 2) {
    Rcpp::stop("biv_iid requires n_pollutant = 2, got %d", n_pol);
  }
  return n_pol;
}

CppAD::AD<double> biv_log_sd(
  const CppAD::vector<CppAD::AD<double>>& x,
  const density_data& model,
  const Config& config,
  int col) {

  const std::size_t t_global = model.theta_index(col);
  const std::size_t t_row = model.theta_row(col);
  if (transform_theta_at(config, t_row)) {
    return x[t_global];
  }
  return CppAD::log(x[t_global]);
}

CppAD::vector<CppAD::AD<double>> biv_iid_dens(
  const CppAD::vector<CppAD::AD<double>>& x,
  const density_data& model,
  const Config& config) {

  const int n_pol = biv_n_pollutant(model);
  const int n_col = model.gamma_map.ncol();
  if (n_col <= 0 || n_col % n_pol != 0) {
    Rcpp::stop(
      "biv_iid gamma columns (%d) must be a positive multiple of %d",
      n_col, n_pol);
  }
  const int n_station = n_col / n_pol;

  const CppAD::AD<double> log_s1 = biv_log_sd(x, model, config, 0);
  const CppAD::AD<double> log_s2 = biv_log_sd(x, model, config, 1);
  const CppAD::AD<double> rho = x[model.theta_index(2)];
  const CppAD::AD<double> s1 = CppAD::exp(log_s1);
  const CppAD::AD<double> s2 = CppAD::exp(log_s2);
  const CppAD::AD<double> one_m_rho2 = CppAD::AD<double>(1.0) - rho * rho;

  CppAD::AD<double> quad = 0.0;
  for (int s = 0; s < n_station; ++s) {
    const std::vector<std::size_t> i1 = model.gamma_global_indices(s * n_pol);
    const std::vector<std::size_t> i2 = model.gamma_global_indices(s * n_pol + 1);
    if (i1.empty() || i2.empty()) {
      Rcpp::stop("biv_iid gamma column has no parameter index");
    }
    const CppAD::AD<double> g1 = x[i1[0]];
    const CppAD::AD<double> g2 = x[i2[0]];
    const CppAD::AD<double> z1 = g1 / s1;
    const CppAD::AD<double> z2 = g2 / s2;
    quad += z1 * z1 + z2 * z2 - CppAD::AD<double>(2.0) * rho * z1 * z2;
  }

  const CppAD::AD<double> log_det_term =
    CppAD::AD<double>(n_station) * (
      log_s1 + log_s2 +
      CppAD::AD<double>(0.5) * CppAD::log(one_m_rho2) +
      CppAD::AD<double>(2.0 * ONEHALFLOGTWOPI));

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -CppAD::AD<double>(0.5) * quad / one_m_rho2 - log_det_term;
  return result;
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> biv_iid_sparsity(
  const density_data& model) {

  const int n_pol = biv_n_pollutant(model);
  const int n_col = model.gamma_map.ncol();
  const int n_station = n_col / n_pol;
  const std::size_t t0 = model.theta_index(0);
  const std::size_t t1 = model.theta_index(1);
  const std::size_t t2 = model.theta_index(2);
  const std::size_t thetas[3] = {t0, t1, t2};

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      pairs.insert({thetas[a], thetas[b]});
    }
  }
  for (int s = 0; s < n_station; ++s) {
    const std::vector<std::size_t> i1 = model.gamma_global_indices(s * n_pol);
    const std::vector<std::size_t> i2 = model.gamma_global_indices(s * n_pol + 1);
    if (i1.empty() || i2.empty()) {
      continue;
    }
    const std::size_t g1 = i1[0];
    const std::size_t g2 = i2[0];
    const std::size_t gs[2] = {g1, g2};
    for (int a = 0; a < 2; ++a) {
      for (int b = 0; b < 2; ++b) {
        pairs.insert({gs[a], gs[b]});
      }
      for (int b = 0; b < 3; ++b) {
        pairs.insert({gs[a], thetas[b]});
        pairs.insert({thetas[b], gs[a]});
      }
    }
  }

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, model.n_tape, pairs);
  return hessian;
}

}  // namespace

//' Build raw AD handle for a bivariate iid station-intercept shard
//'
//' @param model An \code{density_data} object. \code{precision} is a list with
//'   \code{n_pollutant = 2}. Gamma columns are station-major, pollutant-minor.
//'   Theta columns are \code{log(sd1)}, \code{log(sd2)}, and \code{rho}.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP create_ad_shard_biv_iid(SEXP model, Rcpp::List config) {
  const density_data ad_model(model);
  AdTape pack = build_ad_fun_random_with_pattern(
    ad_model, config, biv_iid_dens, biv_iid_sparsity);
  std::vector<AdTape> packs;
  packs.push_back(std::move(pack));
  ad_pack* groups = packs_to_ad_fun(
    std::move(packs), ad_model.num_beta, ad_model.num_theta,
    adlaplace_make_ad_shard);
  return make_ad_pack_ptr(groups);
}
