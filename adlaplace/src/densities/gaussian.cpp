#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/densities/gaussian.hpp"
#include "adlaplace/math/constants.hpp"

// Gaussian observation density y_i ~ N(eta_i, sigma^2) with
// eta = X beta + A gamma and a single scale parameter sigma at theta
// index 0 (log scale when config$transform_theta is TRUE).
//
// The observation shard carries the eta-dependent quadratic term
// -0.5 * (y - eta)^2 / sigma^2 per observation group; the normalizing
// terms -ny * (log(sigma) + 0.5 * log(2*pi)), which depend only on the
// parameters, live in gaussian_extra (see parameter_ad_fun()).
CppAD::vector<CppAD::AD<double>> gaussian_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup) {

  const Config config(config_list);

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> thetaIn = x[theta_index];
  CppAD::AD<double> logSigma = config.transform_theta ? thetaIn : CppAD::log(thetaIn);
  CppAD::AD<double> inv2var = CppAD::AD<double>(0.5) * CppAD::exp(-2 * logSigma);

  const bool have_shards = config.shards.ncol() > 0;
  const size_t ny = model.y.size();
  size_t startP = 0;
  size_t endP = 0;
  if (have_shards) {
    startP = config.shards.p[Dgroup];
    endP = config.shards.p[Dgroup + 1];
  } else if (Dgroup == 0) {
    endP = ny;
  }

  CppAD::AD<double> ss = 0.0;
  for (size_t DI = startP; DI < endP; ++DI) {
    const size_t Dobs = have_shards ? config.shards.i[DI] : DI;

    CppAD::AD<double> eta = 0.0;
    const size_t p0x = model.XTp.p[Dobs];
    const size_t p1x = model.XTp.p[Dobs + 1];
    for (size_t D = p0x; D < p1x; ++D) {
      eta += model.XTp.x[D] * x[model.XTp.i[D]];
    }
    const size_t p0a = model.ATp.p[Dobs];
    const size_t p1a = model.ATp.p[Dobs + 1];
    for (size_t D = p0a; D < p1a; ++D) {
      eta += model.ATp.x[D] * x[model.num_beta + model.ATp.i[D]];
    }

    const CppAD::AD<double> resid = CppAD::AD<double>(model.y[Dobs]) - eta;
    ss += resid * resid;
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -inv2var * ss;
  return result;
}

CppAD::vector<CppAD::AD<double>> gaussian_extra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list) {

  const Config config(config_list);

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> thetaIn = x[theta_index];
  CppAD::AD<double> logSigma = config.transform_theta ? thetaIn : CppAD::log(thetaIn);

  const double ny = static_cast<double>(model.y.size());

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -ny * logSigma - CppAD::AD<double>(ny * ONEHALFLOGTWOPI);
  return result;
}
