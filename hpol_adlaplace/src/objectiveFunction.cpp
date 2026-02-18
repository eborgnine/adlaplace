// #define DEBUG  // DEBUG uses R output and is unsafe with OpenMP worker threads

#include "adlaplace/adlaplace.hpp"
#include "adlaplace/math/lgamma.hpp"


#define COMPUTE_CONSTANTS
//#define DEBUG

// use the standard log density for random effects
//#include "adlaplace/logDens/random.hpp"


CppAD::vector<CppAD::AD<double>> logDensRandom(
  const CppAD::vector<CppAD::AD<double>>& x,
  const Data& data,
  const Config& config
  ){
#ifdef UNDEF
  CppAD::vector<CppAD::AD<double>> logTheta(config.Ntheta), expTheta(config.Ntheta), gammaScaled(config.Ngamma);

  if(config.transform_theta) {
    for(size_t D=0,Dtheta=config.theta_begin;D<config.Ntheta;++D,Dtheta++) {
      logTheta[D] = x[Dtheta];
      expTheta[D] = CppAD::exp(logTheta[D]);
    }
  } else {
    for(size_t D=0,Dtheta=config.theta_begin;D<config.Ntheta;++D,Dtheta++) {
      expTheta[D] = x[Dtheta];
      logTheta[D] = CppAD::log(expTheta[D]);
    }
  }

  CppAD::AD<double> qpart = 0.0, qDet=0.0;

  if(config.verbose) {
    Rcpp::Rcout << "q ngamma  " << data.Ngamma << " map.ncol " << data.map.ncol() << 
    " map@i.size " << data.map.i.size() << 
    " ntheta " << config.theta.size() << 
    " exp theta map 0 " <<  expTheta[data.map.i[0]] << " gamma0 " << x[config.gamma_begin] << ".\n";
  }

  for(size_t D=0,Dgamma=config.gamma_begin;D<config.Ngamma;D++,Dgamma++) {
    gammaScaled[D] = x[Dgamma];
  }

  for(size_t Dtheta=0;Dtheta<data.Nmap;Dtheta++) {
    const size_t mapStart = data.map.p[Dtheta];
    const size_t mapEnd = data.map.p[Dtheta+1]; 
    const size_t Nhere = mapEnd - mapStart;

    if(Nhere) {
      qDet += logTheta[Dtheta] * Nhere;
    }

    for(size_t DgammaI=mapStart;DgammaI<mapEnd;DgammaI++) {
      const size_t Dgamma = data.map.i[DgammaI];
      gammaScaled[Dgamma] /= expTheta[Dtheta];
    }
  }

  for(size_t D=0;D<data.Ngamma;D++) {
    qpart += 
      gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
  }
  qpart *= 0.5;

  // Warning, no offdiag of Q implemented
  qDet += CppAD::AD<double>(data.Ngamma * ONEHALFLOGTWOPI);

#endif  

  CppAD::vector<CppAD::AD<double>> result(1, 0.0);
  result[0] = 0;//- qpart - qDet;

  return(result);
}

#ifdef UNDEF
CppAD::AD<double> stable_logsumexp(const CppAD::vector<CppAD::AD<double>>& eta) {
  // compute log(sum(exp(eta)))
  CppAD::AD<double> max_value_ad = eta[0];
  for (size_t Deta = 1; Deta < eta.size(); ++Deta) {
    max_value_ad = CppAD::CondExpGt(eta[Deta], max_value_ad, eta[Deta], max_value_ad);
  }

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
#endif

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& params,
  const Data& data,
  const Config& config,
  const size_t Dgroup)
{
  CppAD::AD<double> result1 = 0.0;
  CppAD::vector<CppAD::AD<double>> result(1);

#ifdef UNDEF
  const bool have_groups = config.groups.ncol() > 0;
  const size_t startP = have_groups ? config.groups.p[Dgroup] : Dgroup;
  const size_t endP = have_groups ? config.groups.p[Dgroup + 1] : Dgroup + 1;

  const CppAD::AD<double> lastTheta = params[params.size() - 1];
  const CppAD::AD<double> logNuSq = config.transform_theta
    ? 2.0 * lastTheta
    : 2.0 * CppAD::exp(lastTheta);


  for (size_t DstrataI = startP; DstrataI < endP; DstrataI++) {
    const size_t Dstrata = have_groups ? config.groups.i[DstrataI] : DstrataI;

    const auto etaHere = compute_eta_for_stratum(Dstrata, data, config, params);
    const auto contrib = accumulate_contrib_for_stratum(Dstrata, data, etaHere, logNuSq);

    result1 += contrib;
  }
#endif
  result[0] = -result1;
  return result;
}

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& params,
  const Data& data,
  const Config& config)
{
#ifdef UNDEF
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

#endif  
  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = 0;//contrib;

  return result;
}

// ADfun and interfaces
#include "adlaplace/creators/adfun.hpp"
