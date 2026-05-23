#include "adlaplace/densities/random.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/math/constants.hpp"

namespace {

CppAD::vector<CppAD::AD<double>> random_diagonal(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_model& model,
  const NumVecView& Q,
  const std::vector<std::size_t>& gamma_indices,
  const Config& config) {

  const std::size_t Ngamma = gamma_indices.size();
  if (Q.size() != Ngamma) {
    Rcpp::warning("precision$Q length (%d) differs from gamma_map rows (%d)",
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
    Rcpp::Rcout << "theta index " << theta_index <<
      " logVariance " << logSd << " precision " << precision << "\n";
    Rcpp::Rcout << "random_diagonal n_gamma " << Ngamma <<
      " qDet " << qDet <<
      " qpart " << qpart << "\n";
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -qpart - qDet + CppAD::AD<double>(0.5 * logQsum);
  return result;
}

}  // namespace

CppAD::vector<CppAD::AD<double>> random_diagonal(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_model& model,
  const Rcpp::List& precision,
  const Rcpp::List& config) {

  if (!precision.containsElementNamed("Q")) {
    Rcpp::stop("precision must contain element Q for random_diagonal");
  }
  const NumVecView Q(precision["Q"]);
  const std::vector<std::size_t> gamma_indices = model.all_gamma_global_indices();
  return random_diagonal(x, model, Q, gamma_indices, Config(config));
}
