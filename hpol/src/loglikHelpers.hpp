#pragma once

#ifndef HPOL_HELPERS_HPP
#define HPOL_HELPERS_HPP


#include"hpol.hpp"
template <class Type>
Type stable_logsumexp(const CppAD::vector<Type>& eta) {
    // compute log(sum(exp(eta)))

    // find index of maximum
    size_t max_idx = 0;
    for (size_t Deta = 1; Deta < eta.size(); ++Deta) {
        if (eta[Deta] > eta[max_idx]) {
            max_idx = Deta;
        }
    }

    // get max value, sever tape if AD
    double max_value;
    if constexpr (std::is_same<Type, double>::value) {
        max_value = eta[max_idx];
    } else {
        max_value = CppAD::Value(eta[max_idx]);
    }

    // sum exp differences
    Type sumexp = Type(0);
    for (size_t Deta = 0; Deta < eta.size(); ++Deta) {
        sumexp += CppAD::exp(eta[Deta] - max_value);
    }

    return Type(max_value) + CppAD::log(sumexp);
}


// Compute eta for one stratum, using Data directly
template <class GammaT, class BetaT>
CppAD::vector<GammaT> compute_eta_for_stratum(size_t Dstrata,
                             const Data& data,
                             const CppAD::vector<GammaT>& gamma,
                             const CppAD::vector<BetaT>& beta)
{

  const size_t startHere = data.CC.p[Dstrata], endHere = data.CC.p[Dstrata+1];
  const size_t NinStrata = endHere - startHere;

  CppAD::vector<GammaT> etaHere(NinStrata);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const int Deta = data.CC.i[k];

    GammaT accGamma = GammaT(0);
    BetaT accBeta = BetaT(0);

    // X contribution
    const int p0x = data.X.p[Deta];
    const int p1x = data.X.p[Deta + 1];
    for (int t = p0x; t < p1x; ++t) {
      accBeta += data.X.x[t] * beta[data.X.i[t]];
    }

    // A contribution
    const int p0a = data.A.p[Deta];
    const int p1a = data.A.p[Deta + 1];
    for (int t = p0a; t < p1a; ++t) {
      accGamma += data.A.x[t] * gamma[data.A.i[t]];
    }

    etaHere[j] = GammaT(accBeta) + accGamma;
  }
  return(etaHere);
}

// Accumulate contribution over one stratum.
// - Out: accumulation/tape type (double or CppAD::AD<double>)
// - ParamsT: provides logSqrtNu and (optionally) other constants; e.g. PackedParams<double>
// Returns {contrib, sumY}
template <class Out, class ParamsT>
Out accumulate_contrib_for_stratum(size_t Dstrata,
                               const Data& data,
                               const CppAD::vector<Out>& etaHere, // length = NinStrata
                               const PackedParams<ParamsT>& params,
                               const Config& cfg)
{
  using CppAD::exp;

  const size_t startHere = data.CC.p[Dstrata];
  const size_t endHere   = data.CC.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;


  Out   contrib = Out(0);
  int   sumY    = 0;

  if (NinStrata == 0) return contrib;

  // normalize to get probabilities
  Out etaLogSum = stable_logsumexp(etaHere);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const size_t idx = static_cast<size_t>(data.CC.i[k]);

    sumY += static_cast<int>(data.y[idx]);

    const Out etaMinusLogSumMu = etaHere[j] - etaLogSum;
    const Out muBarDivSqrtNu   = exp(etaMinusLogSumMu - Out(params.logSqrtNu));

    if (cfg.dirichlet) {
      contrib += lgamma_ad( Out(data.y[idx]) + muBarDivSqrtNu )
               - lgamma_ad( muBarDivSqrtNu );
    } else {
      contrib += Out(data.y[idx]) * etaMinusLogSumMu;
    }

#ifdef EVALCONSTANTS
    contrib -= lgamma(static_cast<double>(data.y[idx]) + 1.0);
#endif
  }

if (cfg.dirichlet) {
  ParamsT dirichletContrib = params.lgammaOneOverSqrtNu - lgamma_ad(params.oneOverSqrtNu + ParamsT(sumY));
  contrib += Out(dirichletContrib);
}

  return contrib;
}


#endif
