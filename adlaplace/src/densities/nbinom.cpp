#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/math/lgamma.hpp"
#include "adlaplace/densities/nbinom.hpp"

CppAD::vector<CppAD::AD<double>> nbinom_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup) {

  const Config config(config_list);

  CppAD::AD<double> logDens1 = 0.0, logDens2 = 0.0;

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> thetaIn = x[theta_index];
  CppAD::AD<double> logTheta = config.transform_theta ? thetaIn : CppAD::log(thetaIn);

  CppAD::AD<double> logNbSize = -2 * logTheta;
  CppAD::AD<double> nbSize = CppAD::exp(logNbSize);

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

  for (size_t DI = startP; DI < endP; DI++) {
    const size_t Dobs = have_shards ? config.shards.i[DI] : DI;

    const size_t p0x = model.XTp.p[Dobs];
    const size_t p1x = model.XTp.p[Dobs + 1];
    const size_t p0a = model.ATp.p[Dobs];
    const size_t p1a = model.ATp.p[Dobs + 1];
    CppAD::AD<double> etaFixed = 0.0;
    for (size_t D = p0x; D < p1x; D++) {
      etaFixed += model.XTp.x[D] * x[model.XTp.i[D]];
    }
    CppAD::AD<double> etaRandom = 0.0;
    for (size_t D = p0a; D < p1a; D++) {
      etaRandom += model.ATp.x[D] * x[model.num_beta + model.ATp.i[D]];
    }

    const CppAD::AD<double> eta = etaRandom + etaFixed;

    const CppAD::AD<double> diff = eta - logNbSize;
    const CppAD::AD<double> softplus = CppAD::CondExpGt(
      diff, CppAD::AD<double>(0.0),
      diff + CppAD::log1p(CppAD::exp(-diff)),
      CppAD::log1p(CppAD::exp(diff)));

    const CppAD::AD<double> logRplusMu = logNbSize + softplus;

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
  const Rcpp::List& config_list) {

  const Config config(config_list);

  const std::size_t theta_index = model.theta_index(0);
  CppAD::AD<double> thetaIn = x[theta_index];
  CppAD::AD<double> logTheta = config.transform_theta ? thetaIn : CppAD::log(thetaIn);

  CppAD::AD<double> logNbSize = -2 * logTheta;
  CppAD::AD<double> nbSize = CppAD::exp(logNbSize);
  CppAD::AD<double> lgammaNbSize = lgamma_ad(nbSize);
  CppAD::AD<double> sizeLogSize = nbSize * logNbSize;

  const std::size_t ny = model.y.size();
  CppAD::AD<double> logDens1 = 0.0;

  for (std::size_t D = 0; D < ny; D++) {
    logDens1 += lgamma_ad(model.y[D] + nbSize);
  }

  CppAD::AD<double> logDens2 = static_cast<double>(ny) * (sizeLogSize - lgammaNbSize);

  double constants = 0.0;
  for (std::size_t D = 0; D < ny; D++) {
    constants += std::lgamma(model.y[D] + 1.0);
  }
  logDens2 -= constants;

  CppAD::AD<double> logDens = logDens1 + logDens2;
  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = logDens;
  return result;
}
