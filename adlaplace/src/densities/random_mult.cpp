#include "adlaplace/densities/random_mult.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/math/constants.hpp"

CppAD::vector<CppAD::AD<double>>
random_mult(const CppAD::vector<CppAD::AD<double>> &x, const ad_data &model,
            SEXP precision, const Rcpp::List &config_list) {

  const Config config(config_list);

  if (Rf_isNull(precision) || TYPEOF(precision) != VECSXP) {
    Rcpp::stop(
        "random_mult precision must be list(Q = <dgCMatrix>, log_det, rank)");
  }
  Rcpp::List prec(precision);
  if (!prec.containsElementNamed("Q")) {
    Rcpp::stop("random_mult precision list must contain Q");
  }
  SEXP q_sexp = prec["Q"];
  if (!Rf_inherits(q_sexp, "dgCMatrix") && !Rf_inherits(q_sexp, "ngCMatrix")) {
    Rcpp::stop("random_mult precision Q must be a dgCMatrix "
               "(convert with as(Q, \"generalMatrix\") for symmetric storage)");
  }
  const DgCView Q(Rcpp::as<Rcpp::S4>(q_sexp));

  const int n_term = model.gamma_map.ncol();
  if (static_cast<std::size_t>(Q.nrow()) != static_cast<std::size_t>(n_term) ||
      static_cast<std::size_t>(Q.ncol()) != static_cast<std::size_t>(n_term)) {
    Rcpp::stop("random_mult Q is %d x %d but gamma_map has %d columns",
               Q.nrow(), Q.ncol(), n_term);
  }

  const double rank = prec.containsElementNamed("rank")
                          ? Rcpp::as<double>(prec["rank"])
                          : static_cast<double>(n_term);
  const double log_det = prec.containsElementNamed("log_det")
                             ? Rcpp::as<double>(prec["log_det"])
                             : 0.0;

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> logSd;
  if (config.transform_theta) {
    logSd = x[theta_index];
  } else {
    logSd = CppAD::log(x[theta_index]);
  }
  CppAD::AD<double> tau = CppAD::exp(-2 * logSd);

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
  CppAD::AD<double> qpart = CppAD::AD<double>(0.5) * tau * qf;

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
