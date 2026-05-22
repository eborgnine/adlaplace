#include "adlaplace/densities/random.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/math/constants.hpp"

namespace {

struct data_random_diag {
  NumVecView Q;
  IntVecView gamma_map;
  size_t theta_index;
  size_t Ngamma;

  explicit data_random_diag(const Rcpp::List& data)
    : Q(data["Q"]),
      gamma_map(data["gamma_map"]) {
    Ngamma = gamma_map.size();
    IntVecView theta_map(data["theta_map"]);
    if (theta_map.size() > 0) {
      theta_index = static_cast<size_t>(theta_map[0]);
    }
    if (gamma_map.size() != Q.size()) {
      Rcpp::warning("gamma_map and Q have different lengths in data_random_diag");
    }
  }
};

CppAD::vector<CppAD::AD<double>> random_diagonal(
  const CppAD::vector<CppAD::AD<double>>& x,
  const data_random_diag& data,
  const Config& config) {

  CppAD::AD<double> logSd;

  // theta is SD
  if (config.transform_theta) {
    logSd = x[data.theta_index];
  } else {
    logSd = CppAD::log(x[data.theta_index]);
  }

  CppAD::AD<double> precision = CppAD::exp(-2 * logSd);
  CppAD::AD<double> qpart = 0.0;
  double logQsum = 0.0;

  for (size_t k = 0; k < data.Ngamma; ++k) {
    const size_t gidx = static_cast<size_t>(data.gamma_map[k]);
    const double qk = data.Q[k];
    qpart += x[gidx] * x[gidx] * qk;
    logQsum += std::log(qk);
  }
  qpart *= CppAD::AD<double>(0.5) * precision;

  CppAD::AD<double> qDet = logSd * CppAD::AD<double>(data.Ngamma) +
    CppAD::AD<double>(data.Ngamma * ONEHALFLOGTWOPI);

    if(config.verbose) {
      Rcpp::Rcout << "theta index " << data.theta_index << 
      " logVariance " << logSd << " precision " << precision << "\n";
      Rcpp::Rcout << "random_diagonal n_gamma " << data.Ngamma << 
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
  const Rcpp::List& data,
  const Rcpp::List& config) {
  return random_diagonal(x, data_random_diag(data), Config(config));
}
