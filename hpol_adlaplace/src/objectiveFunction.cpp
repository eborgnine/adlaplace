#define DEBUG

#include "adlaplace/adlaplace.hpp"
//#include "adlaplace/math/lgamma.hpp"

CppAD::AD<double> lgamma_ad(const CppAD::AD<double>& x) {
  return x;
}

#define COMPUTE_CONSTANTS
//#define DEBUG

// use the standard log density for random effects
#include "adlaplace/logDens/random.hpp"

CppAD::AD<double> stable_logsumexp(const CppAD::vector<CppAD::AD<double>>& eta) {
  // compute log(sum(exp(eta)))

  size_t max_idx = 0;
  for (size_t Deta = 1; Deta < eta.size(); ++Deta) {
    if (eta[Deta] > eta[max_idx]) {
      max_idx = Deta;
    }
  }

  // sever tape for the max shift
  const double max_value = CppAD::Value(eta[max_idx]);
  const CppAD::AD<double> max_value_ad = max_value;

  CppAD::AD<double> sumexp = 0.0;
  for (size_t Deta = 0; Deta < eta.size(); ++Deta) {
    sumexp += CppAD::exp(eta[Deta] - max_value_ad);
  }

  return CppAD::log(sumexp) + max_value_ad;
}

// Compute eta for one stratum using the full parameter vector.
CppAD::vector<CppAD::AD<double>> compute_eta_for_stratum(
  const size_t Dstrata,
  const Data& data,
  const Config& config,
  const CppAD::vector<CppAD::AD<double>>& params)
{
  const size_t startHere = data.elgm_matrix.p[Dstrata];
  const size_t endHere = data.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  CppAD::vector<CppAD::AD<double>> etaHere(NinStrata);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const int Deta = data.elgm_matrix.i[k];

    CppAD::AD<double> accGamma = 0.0;
    CppAD::AD<double> accBeta = 0.0;

    // X contribution
    const int p0x = data.X.p[Deta];
    const int p1x = data.X.p[Deta + 1];
    for (int t = p0x; t < p1x; ++t) {
      accBeta += data.X.x[t] * params[config.beta_begin + data.X.i[t]];
    }

    // A contribution
    const int p0a = data.A.p[Deta];
    const int p1a = data.A.p[Deta + 1];
    for (int t = p0a; t < p1a; ++t) {
      accGamma += data.A.x[t] * params[config.gamma_begin + data.A.i[t]];
    }

    etaHere[j] = accBeta + accGamma;
  }

  return etaHere;
}

// Accumulate contribution over one stratum.
CppAD::AD<double> accumulate_contrib_for_stratum(
  const size_t Dstrata,
  const Data& data,
  const CppAD::vector<CppAD::AD<double>>& etaHere,
  const CppAD::AD<double>& logNuSq)
{
  const size_t startHere = data.elgm_matrix.p[Dstrata];
  const size_t endHere = data.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  CppAD::AD<double> contrib = 0.0;
  if (NinStrata == 0) {
    return contrib;
  }

  // normalize to get probabilities
  const CppAD::AD<double> etaLogSum = stable_logsumexp(etaHere);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const size_t idx = static_cast<size_t>(data.elgm_matrix.i[k]);
    const double yhere = data.y[idx];

    const CppAD::AD<double> etaMinusLogSumMu = etaHere[j] - etaLogSum;
    const CppAD::AD<double> muBarDivNuSq = CppAD::exp(etaMinusLogSumMu - logNuSq);

    contrib += lgamma_ad(muBarDivNuSq + yhere) - lgamma_ad(muBarDivNuSq);
  }

  return contrib;
}

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& params,
  const Data& data,
  const Config& config,
  const size_t Dgroup)
{
  const bool have_groups = config.groups.ncol() > 0;
  const size_t startP = have_groups ? config.groups.p[Dgroup] : Dgroup;
  const size_t endP = have_groups ? config.groups.p[Dgroup + 1] : Dgroup + 1;

  const CppAD::AD<double> lastTheta = params[params.size() - 1];
  const CppAD::AD<double> logNuSq = config.transform_theta
    ? 2.0 * lastTheta
    : 2.0 * CppAD::exp(lastTheta);

  CppAD::AD<double> result1 = 0.0;
  CppAD::vector<CppAD::AD<double>> result(1);

  for (size_t DstrataI = startP; DstrataI < endP; DstrataI++) {
    const size_t Dstrata = have_groups ? config.groups.i[DstrataI] : DstrataI;

    const auto etaHere = compute_eta_for_stratum(Dstrata, data, config, params);
    const auto contrib = accumulate_contrib_for_stratum(Dstrata, data, etaHere, logNuSq);

    result1 += contrib;
  }

  result[0] = -result1;
  return result;
}

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& params,
  const Data& data,
  const Config& config)
{
  const CppAD::AD<double> lastTheta = params[params.size() - 1];
  const CppAD::AD<double> logNuSq = config.transform_theta
    ? 2.0 * lastTheta
    : 2.0 * CppAD::exp(lastTheta);
  const CppAD::AD<double> oneOverNuSq = CppAD::exp(-logNuSq);

  const size_t Nstrata = data.elgm_matrix.ncol();

  CppAD::AD<double> contribLgamma1overNuSqPlusSumYil = 0.0;
#ifdef COMPUTE_CONSTANTS
  CppAD::AD<double> contribLgamma1pSumYil = 0.0;
  CppAD::AD<double> contribLgammaYp1 = 0.0;
#endif

  for (size_t Dstrata = 0; Dstrata < Nstrata; Dstrata++) {
    int sumYhere = 0;
    const size_t endHere = data.elgm_matrix.p[Dstrata + 1];

    for (size_t DobsI = data.elgm_matrix.p[Dstrata]; DobsI < endHere; DobsI++) {
      const size_t Dobs = data.elgm_matrix.i[DobsI];
      sumYhere += data.y[Dobs];
#ifdef COMPUTE_CONSTANTS
      contribLgammaYp1 += std::lgamma(1 + data.y[Dobs]);
#endif
    }

    contribLgamma1overNuSqPlusSumYil += lgamma_ad(oneOverNuSq + sumYhere);
#ifdef COMPUTE_CONSTANTS
    contribLgamma1pSumYil += std::lgamma(1 + sumYhere);
#endif
  }

  CppAD::AD<double> contrib =
    Nstrata * lgamma_ad(oneOverNuSq) + contribLgamma1overNuSqPlusSumYil;

#ifdef COMPUTE_CONSTANTS
  contrib += contribLgammaYp1 + contribLgamma1pSumYil;
#endif
  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = contrib;
  return result;
}

// ADfun and interfaces
#include "adlaplace/creators/adfun.hpp"
