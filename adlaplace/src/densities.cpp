#include "adlaplace/ad_data.hpp"
#include "adlaplace/atomics.hpp"
#include "adlaplace/densities.hpp"
#include "adlaplace/eta.hpp"
#include "adlaplace/rviews.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace {

std::pair<std::size_t, std::size_t> obs_range(
  const Config& config,
  std::size_t ny,
  std::size_t Dgroup) {

  const bool have_shards = config.shards.ncol() > 0;
  std::size_t startP = 0;
  std::size_t endP = 0;
  if (have_shards) {
    startP = static_cast<std::size_t>(config.shards.p[Dgroup]);
    endP = static_cast<std::size_t>(config.shards.p[Dgroup + 1]);
  } else if (Dgroup == 0) {
    endP = ny;
  }
  return {startP, endP};
}

CppAD::AD<double> softplus_ad(const CppAD::AD<double>& z) {
  return CppAD::CondExpGt(
    z, CppAD::AD<double>(0.0),
    z + CppAD::log1p(CppAD::exp(-z)),
    CppAD::log1p(CppAD::exp(z)));
}

CppAD::AD<double> log_theta_ad(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config) {

  const std::size_t theta_index = model.theta_index(0);
  const CppAD::AD<double> thetaIn = x[theta_index];
  return config.transform_theta ? thetaIn : CppAD::log(thetaIn);
}

}  // namespace

CppAD::vector<CppAD::AD<double>> binomial_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config,
  const size_t Dgroup) {

  const size_t ny = model.y.size();
  const auto range = obs_range(config, ny, Dgroup);
  const bool have_shards = config.shards.ncol() > 0;

  CppAD::AD<double> result = 0.0;
  for (size_t DI = range.first; DI < range.second; ++DI) {
    const size_t Dobs = have_shards ? static_cast<size_t>(config.shards.i[DI]) : DI;
    const CppAD::AD<double> eta = eta_at(x, model, Dobs);
    result += CppAD::AD<double>(model.y[Dobs]) * eta - softplus_ad(eta);
  }

  CppAD::vector<CppAD::AD<double>> out(1);
  out[0] = result;
  return out;
}

CppAD::vector<CppAD::AD<double>> gaussian_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config,
  const size_t Dgroup) {

  const CppAD::AD<double> logSigma = log_theta_ad(x, model, config);
  const CppAD::AD<double> inv2var = CppAD::AD<double>(0.5) * CppAD::exp(-2 * logSigma);

  const size_t ny = model.y.size();
  const auto range = obs_range(config, ny, Dgroup);
  const bool have_shards = config.shards.ncol() > 0;

  CppAD::AD<double> ss = 0.0;
  for (size_t DI = range.first; DI < range.second; ++DI) {
    const size_t Dobs = have_shards ? static_cast<size_t>(config.shards.i[DI]) : DI;
    const CppAD::AD<double> eta = eta_at(x, model, Dobs);
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
  const Config& config) {

  const CppAD::AD<double> logSigma = log_theta_ad(x, model, config);
  const double ny = static_cast<double>(model.y.size());

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -ny * logSigma - CppAD::AD<double>(ny * ONEHALFLOGTWOPI);
  return result;
}

CppAD::vector<CppAD::AD<double>> nbinom_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config,
  const size_t Dgroup) {

  CppAD::AD<double> logDens1 = 0.0;
  CppAD::AD<double> logDens2 = 0.0;

  const CppAD::AD<double> logTheta = log_theta_ad(x, model, config);
  const CppAD::AD<double> logNbSize = -2 * logTheta;
  const CppAD::AD<double> nbSize = CppAD::exp(logNbSize);

  const size_t ny = model.y.size();
  const auto range = obs_range(config, ny, Dgroup);
  const bool have_shards = config.shards.ncol() > 0;

  for (size_t DI = range.first; DI < range.second; ++DI) {
    const size_t Dobs = have_shards ? static_cast<size_t>(config.shards.i[DI]) : DI;
    const CppAD::AD<double> eta = eta_at(x, model, Dobs);
    const CppAD::AD<double> diff = eta - logNbSize;
    const CppAD::AD<double> logRplusMu = logNbSize + softplus_ad(diff);

    logDens1 += model.y[Dobs] * eta;
    logDens2 += -logRplusMu * (model.y[Dobs] + nbSize);
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = logDens1 + logDens2;
  return result;
}

CppAD::vector<CppAD::AD<double>> nbinom_extra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config) {

  const CppAD::AD<double> logTheta = log_theta_ad(x, model, config);
  const CppAD::AD<double> logNbSize = -2 * logTheta;
  const CppAD::AD<double> nbSize = CppAD::exp(logNbSize);
  const CppAD::AD<double> lgammaNbSize = lgamma_ad(nbSize);
  const CppAD::AD<double> sizeLogSize = nbSize * logNbSize;

  const std::size_t ny = model.y.size();
  CppAD::AD<double> logDens1 = 0.0;
  for (std::size_t D = 0; D < ny; ++D) {
    logDens1 += lgamma_ad(model.y[D] + nbSize);
  }

  CppAD::AD<double> logDens2 = static_cast<double>(ny) * (sizeLogSize - lgammaNbSize);

  double constants = 0.0;
  for (std::size_t D = 0; D < ny; ++D) {
    constants += std::lgamma(model.y[D] + 1.0);
  }
  logDens2 -= constants;

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = logDens1 + logDens2;
  return result;
}

CppAD::vector<CppAD::AD<double>> poisson_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config,
  const size_t Dgroup) {

  const bool have_offset = !config.offset.empty();
  if (have_offset &&
      config.offset.size() != 1 &&
      config.offset.size() != model.y.size()) {
    Rcpp::stop(
      "config$offset must have length 1 or length(y) (%d), got %d",
      static_cast<int>(model.y.size()), static_cast<int>(config.offset.size()));
  }

  const size_t ny = model.y.size();
  const auto range = obs_range(config, ny, Dgroup);
  const bool have_shards = config.shards.ncol() > 0;

  CppAD::AD<double> linear = 0.0;
  CppAD::AD<double> expsum = 0.0;
  double constants = 0.0;

  for (size_t DI = range.first; DI < range.second; ++DI) {
    const size_t Dobs = have_shards ? static_cast<size_t>(config.shards.i[DI]) : DI;
    const CppAD::AD<double> eta = eta_at(x, model, Dobs);

    const double off = have_offset
      ? (config.offset.size() == 1 ? config.offset[0] : config.offset[Dobs])
      : 0.0;
    const double y = model.y[Dobs];

    linear += y * eta;
    expsum += CppAD::exp(eta + off);
    constants += y * off - std::lgamma(y + 1.0);
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = linear - expsum + CppAD::AD<double>(constants);
  return result;
}
