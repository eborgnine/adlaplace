// #define DEBUG  // DEBUG uses R output and is unsafe with OpenMP worker threads

#include "adlaplace/adlaplace.hpp"

// use the standard log density for random effects
#include "adlaplace/logDens/random.hpp"

CppAD::AD<double> stable_logsumexp(const CppAD::vector<CppAD::AD<double>> &eta)
{
  // compute log(sum(exp(eta)))
#ifdef UNDEF
  CppAD::AD<double> max_value_ad = eta[0];
  for (size_t Deta = 1; Deta < eta.size(); ++Deta)
  {
    max_value_ad = CppAD::CondExpGt(eta[Deta], max_value_ad, eta[Deta], max_value_ad);
  }
#endif
  // Find the maximum value (not part of the AD tape)
  size_t max_value = CppAD::Value(eta[0]);
  for (size_t Deta = 1; Deta < eta.size(); ++Deta) {
    double current_value = CppAD::Value(eta[Deta]);
    if (current_value > max_value) {
      max_value = current_value;
    }
  }

  CppAD::AD<double> sumexp = 0.0;
  for (size_t Deta = 0; Deta < eta.size(); ++Deta)
  {
    sumexp += CppAD::exp(eta[Deta] - max_value);
  }

  return CppAD::log(sumexp) + max_value;
}

// Compute eta for one stratum using the full parameter vector.
CppAD::vector<CppAD::AD<double>> compute_eta_for_stratum(
    const size_t Dstrata,
    const Data &data,
    const Config &config,
    const CppAD::vector<CppAD::AD<double>> &params)
{
  const size_t startHere = data.elgm_matrix.p[Dstrata];
  const size_t endHere = data.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  CppAD::vector<CppAD::AD<double>> etaHere(NinStrata);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k)
  {
    const int Deta = data.elgm_matrix.i[k];

    CppAD::AD<double> accGamma = 0.0;
    CppAD::AD<double> accBeta = 0.0;

    // X contribution
    const int p0x = data.X.p[Deta];
    const int p1x = data.X.p[Deta + 1];
    for (int t = p0x; t < p1x; ++t)
    {
      accBeta += data.X.x[t] * params[config.beta_begin + data.X.i[t]];
    }

    // A contribution
    const int p0a = data.A.p[Deta];
    const int p1a = data.A.p[Deta + 1];
    for (int t = p0a; t < p1a; ++t)
    {
      accGamma += data.A.x[t] * params[config.gamma_begin + data.A.i[t]];
    }

    etaHere[j] = accBeta + accGamma;
  }

  return etaHere;
}


CppAD::AD<double> accumulate_contrib_for_stratum(
    const size_t Dstrata,
    const Data &data,
    const CppAD::vector<CppAD::AD<double>> &etaHere,
    const CppAD::AD<double> &tauSq)
{
  const size_t startHere = data.elgm_matrix.p[Dstrata];
  const size_t endHere = data.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  // normalize to get probabilities, etaLogSum = sum(exp(eta))
  const CppAD::AD<double> etaLogSum = stable_logsumexp(etaHere);

  CppAD::AD<double> contrib = 0.0;
  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k)
  {
    const size_t idx = static_cast<size_t>(data.elgm_matrix.i[k]);
    const double yhere = data.y[idx];

    const CppAD::AD<double> logMuKhere = etaHere[j] - etaLogSum;

    // \sum_{i=0}^{y_k-1} \log(\mu_k + i \tau^2)
    if(yhere) {
      contrib += logMuKhere;
    }
    if(yhere > 1) {
      const CppAD::AD<double> expMuKhere = CppAD::exp(logMuKhere);
      for (int i = 1; i < yhere; ++i) {
        contrib += CppAD::log(expMuKhere + i*tauSq);
      }
    }
  } // loop through obs in strata

  return contrib;
}

CppAD::vector<CppAD::AD<double>> logDensObs(
    const CppAD::vector<CppAD::AD<double>> &params,
    const Data &data,
    const Config &config,
    const size_t Dgroup)
{
  CppAD::AD<double> result1 = 0.0;
  CppAD::vector<CppAD::AD<double>> result(1);

  const bool have_groups = config.groups.ncol() > 0;
  const size_t startP = have_groups ? config.groups.p[Dgroup] : Dgroup;
  const size_t endP = have_groups ? config.groups.p[Dgroup + 1] : Dgroup + 1;

  const CppAD::AD<double> lastTheta = params[params.size() - 1]; // log(\tau) if transform_theta otherwise tau
  const CppAD::AD<double> tauSq = config.transform_theta
                                        ? CppAD::exp(2.0 * lastTheta)
                                        : lastTheta * lastTheta;

  for (size_t DstrataI = startP; DstrataI < endP; DstrataI++)
  {
    const size_t Dstrata = have_groups ? config.groups.i[DstrataI] : DstrataI;

    const auto etaHere = compute_eta_for_stratum(Dstrata, data, config, params);
    const auto contrib = accumulate_contrib_for_stratum(Dstrata, data, etaHere, tauSq);

    result1 += contrib;
  }
  result[0] = result1;
  return result;
}

CppAD::vector<CppAD::AD<double>> logDensExtra(
    const CppAD::vector<CppAD::AD<double>> &params,
    const Data &data,
    const Config &config)
{

// \log(N!) - \sum_{k=1}^K \log(y_k!) - \sum_{j=0}^{N-1} \log(1 + j \tau^2) 

  const CppAD::AD<double> lastTheta = params[params.size() - 1];
  const CppAD::AD<double> tauSq = config.transform_theta
                                        ? CppAD::exp(2.0 * lastTheta)
                                        : lastTheta * lastTheta;

  const size_t Nstrata = data.elgm_matrix.ncol();

  CppAD::AD<double> contribLog1pjTauSq = 0.0;
  double contribLogYkFact = 0.0;
  double contribLogNfact = 0.0;

  for (size_t Dstrata = 0; Dstrata < Nstrata; Dstrata++)
  {
    int sumYhere = 0;
    const size_t endHere = data.elgm_matrix.p[Dstrata + 1];

    for (size_t DobsI = data.elgm_matrix.p[Dstrata]; DobsI < endHere; DobsI++)
    {
      const size_t Dobs = data.elgm_matrix.i[DobsI];
      sumYhere += data.y[Dobs];
      contribLogYkFact += std::lgamma(1 + data.y[Dobs]);
    }

    contribLogNfact += std::lgamma(1 + sumYhere); 
    for(size_t Dsum=1;Dsum < sumYhere;Dsum++) {
      contribLog1pjTauSq += CppAD::log1p(Dsum * tauSq);
    }

  }

  CppAD::AD<double> contrib = contribLogNfact - contribLogYkFact - contribLog1pjTauSq;

  if (config.verbose)
  {
    Rcpp::Rcout << "logDensExtra " << contrib << " tauSq " << tauSq << " contribLogNfact " << 
      contribLogNfact << " contribLogYkFact " << contribLogYkFact << " contribLog1pjTauSq " << contribLog1pjTauSq << 
      "\n";
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = contrib;

  return result;
}


#ifdef USE_OLD_STUFF

#include "adlaplace/math/lgamma.hpp"

CppAD::AD<double> accumulate_contrib_for_stratum_old(
    const size_t Dstrata,
    const Data &data,
    const CppAD::vector<CppAD::AD<double>> &etaHere,
    const CppAD::AD<double> &logNuSq)
{
  const size_t startHere = data.elgm_matrix.p[Dstrata];
  const size_t endHere = data.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  CppAD::AD<double> contrib = 0.0;
  if (NinStrata == 0)
  {
    return contrib;
  }

  // normalize to get probabilities
  const CppAD::AD<double> etaLogSum = stable_logsumexp(etaHere);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k)
  {
    const size_t idx = static_cast<size_t>(data.elgm_matrix.i[k]);
    const double yhere = data.y[idx];

    const CppAD::AD<double> etaMinusLogSumMu = etaHere[j] - etaLogSum;
    const CppAD::AD<double> muBarDivNuSq = CppAD::exp(etaMinusLogSumMu - logNuSq);

    contrib += lgamma_ad(muBarDivNuSq + yhere) - lgamma_ad(muBarDivNuSq);
  }

  return contrib;
}


CppAD::vector<CppAD::AD<double>> logDensObs_old(
    const CppAD::vector<CppAD::AD<double>> &params,
    const Data &data,
    const Config &config,
    const size_t Dgroup)
{
  CppAD::AD<double> result1 = 0.0;
  CppAD::vector<CppAD::AD<double>> result(1);

  const bool have_groups = config.groups.ncol() > 0;
  const size_t startP = have_groups ? config.groups.p[Dgroup] : Dgroup;
  const size_t endP = have_groups ? config.groups.p[Dgroup + 1] : Dgroup + 1;

  const CppAD::AD<double> lastTheta = params[params.size() - 1];
  const CppAD::AD<double> logNuSq = config.transform_theta
                                        ? 2.0 * lastTheta
                                        : 2.0 * CppAD::exp(lastTheta);

  for (size_t DstrataI = startP; DstrataI < endP; DstrataI++)
  {
    const size_t Dstrata = have_groups ? config.groups.i[DstrataI] : DstrataI;

    const auto etaHere = compute_eta_for_stratum(Dstrata, data, config, params);
    const auto contrib = accumulate_contrib_for_stratum_old(Dstrata, data, etaHere, logNuSq);

    result1 += contrib;
  }
  result[0] = result1;
  return result;
}

CppAD::vector<CppAD::AD<double>> logDensExtra_old(
    const CppAD::vector<CppAD::AD<double>> &params,
    const Data &data,
    const Config &config)
{
  // nu is the SD
  // 1/nu^2 is the gamma shape (or precision)

  const CppAD::AD<double> lastTheta = params[params.size() - 1];
  const CppAD::AD<double> logNuSq = config.transform_theta
                                        ? 2.0 * lastTheta
                                        : 2.0 * CppAD::exp(lastTheta);
  const CppAD::AD<double> oneOverNuSq = CppAD::exp(-logNuSq);

  const size_t Nstrata = data.elgm_matrix.ncol();

  CppAD::AD<double> contribLgamma1overNuSqPlusSumYil = 0.0;
  double contribLgamma1pSumYil = 0.0;
  double contribLgammaYp1 = 0.0;
  CppAD::AD<double> lgammaOneOverNuSq = lgamma_ad(oneOverNuSq);

  for (size_t Dstrata = 0; Dstrata < Nstrata; Dstrata++)
  {
    int sumYhere = 0;
    const size_t endHere = data.elgm_matrix.p[Dstrata + 1];

    for (size_t DobsI = data.elgm_matrix.p[Dstrata]; DobsI < endHere; DobsI++)
    {
      const size_t Dobs = data.elgm_matrix.i[DobsI];
      sumYhere += data.y[Dobs];
      contribLgammaYp1 += std::lgamma(1 + data.y[Dobs]);
    }

    contribLgamma1overNuSqPlusSumYil += lgamma_ad(oneOverNuSq + sumYhere);
    contribLgamma1pSumYil += std::lgamma(1 + sumYhere);
  }

  CppAD::AD<double> contribSumN = 0.0;
  for (size_t Dstrata = 0; Dstrata < Nstrata; Dstrata++)
  {
    int sumYhere = 0;
    const size_t endHere = data.elgm_matrix.p[Dstrata + 1];

    for (size_t DobsI = data.elgm_matrix.p[Dstrata]; DobsI < endHere; DobsI++)
    {
      const size_t Dobs = data.elgm_matrix.i[DobsI];
      sumYhere += data.y[Dobs];
    }

    for (int j = 0; j < sumYhere; ++j) {
      contribSumN -= CppAD::log(1.0 + j * CppAD::exp(-logNuSq));
    }
  }

  CppAD::AD<double> contrib =
      Nstrata * lgammaOneOverNuSq - contribLgamma1overNuSqPlusSumYil - contribLgammaYp1 + contribLgamma1pSumYil + contribSumN;

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = contrib;

  if (config.verbose)
  {
    Rcpp::Rcout << "logDensExtra " << contrib << " lasttheta " << lastTheta << " oneOverNuSq " << oneOverNuSq << " logNuSq " << logNuSq << " lgammaOneOverNuSq " << lgammaOneOverNuSq << " contribLgamma1overNuSqPlusSumYil " << contribLgamma1overNuSqPlusSumYil << " contribLgamma1pSumYil " << contribLgamma1pSumYil << " contribLgammaYp1 " << contribLgammaYp1 << "\n";
  }

  return result;
}
#endif

// ADfun and interfaces
#include "adlaplace/creators/adfun.hpp"



