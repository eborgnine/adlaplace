#pragma once

#ifndef HPOL_HELPERS_HPP
#define HPOL_HELPERS_HPP

#include"lgamma.hpp"

//#define DEBUG

//#include"hpol.hpp"
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
    Type max_valueT = Type(max_value);

    for (size_t Deta = 0; Deta < eta.size(); ++Deta) {
        sumexp += CppAD::exp(eta[Deta] - max_valueT);
    }

    Type logSum = CppAD::log(sumexp);
    Type result = logSum + max_valueT;

    return(result);
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
      accGamma += GammaT(data.A.x[t]) * gamma[data.A.i[t]];
    }

    etaHere[j] = GammaT(accBeta) + accGamma;
  }
  return(etaHere);
}

// Accumulate contribution over one stratum.
// - Out: accumulation/tape type (double or CppAD::AD<double>)
// - ParamsT: provides logNuSq and (optionally) other constants; e.g. PackedParams<double>
// Returns {contrib, sumY}
template <class Out, class ParamsT>
CppAD::vector<Out> accumulate_contrib_for_stratum(size_t Dstrata,
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
  CppAD::vector<Out> result(1);

  if (NinStrata == 0) {
    result[0] = contrib;
    return(result);
  }

  // normalize to get probabilities
  Out etaLogSum = stable_logsumexp<Out>(etaHere);


  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const size_t idx = static_cast<size_t>(data.CC.i[k]);
    const auto yhere = data.y[idx];
    Out yOut = Out(yhere);

    sumY += yhere;

    const Out etaMinusLogSumMu = etaHere[j] - etaLogSum;
    const Out muBarDivNuSq   = cfg.dirichlet ?  CppAD::exp(etaMinusLogSumMu - Out(params.logNuSq)):0;

    if (cfg.dirichlet) {
      contrib += lgamma_any<Out>(yOut + muBarDivNuSq) - lgamma_any<Out>(muBarDivNuSq);
      contrib -= Out(lgamma_any<Out>(static_cast<double>(data.y[idx]) + 1.0));
#ifdef DEBUG
      Rcpp::Rcout << lgamma_any<Out>(yOut + muBarDivNuSq)  << " m " << muBarDivNuSq << " " << 
        lgamma_any<Out>(muBarDivNuSq) << " " <<  data.y[idx] + 1.0 << " " << lgamma_any<Out>(static_cast<double>(data.y[idx]) + 1.0) <<
        "\n";
#endif      
    } else {
      contrib += yOut * etaMinusLogSumMu;
    }
  }

if (cfg.dirichlet) {
  ParamsT dirichletContrib = ParamsT(params.lgammaOneOverNuSq) - 
    lgamma_any<ParamsT>(params.oneOverNuSq + sumY) +
    lgamma_any<ParamsT>(1 + sumY);
#ifdef DEBUG
      Rcpp::Rcout << "dirichletContrib " << dirichletContrib <<
      " " <<  ParamsT(params.lgammaOneOverNuSq)  << " " <<  lgamma_any<ParamsT>(params.oneOverNuSq + sumY) 
      << " " << lgamma_any<ParamsT>(1 + sumY) << "\n";
#endif      

  contrib += Out(dirichletContrib);
}

  result[0] = contrib;
  return result;
}


template<class Type, class ParamsT> CppAD::vector<Type> loglikOneStrata(
const int Dstrata,
const CppAD::vector<Type>& gamma,
const PackedParams<ParamsT>& parameters, // gamma is ignored
const Data& data, 
const Config& cfg
){


  auto etaHere = compute_eta_for_stratum<Type, ParamsT>(
    Dstrata, data, gamma, parameters.beta);

  auto contrib = accumulate_contrib_for_stratum<Type, ParamsT>(
    Dstrata, data, etaHere, parameters, cfg
    );

  return(contrib);
}

template<class Type, class ParamsT> CppAD::vector<Type> loglikNStrata(
const int Dstrata,
const int Niter,
const CppAD::vector<Type>& gamma,
const PackedParams<ParamsT>& parameters, // gamma is ignored
const Data& data, 
const Config& cfg
){
  const size_t end = data.Nstrata < (Dstrata+Niter) ? data.Nstrata : (Dstrata+Niter);
  CppAD::vector<Type> result(1);
  result[0] = Type(0);

  for(size_t Diter=Dstrata; Diter < end; ++Diter ) {
    result[0] += loglikOneStrata(Diter, gamma, parameters, data, cfg)[0];
  }
  return(result);
}

template<class Type, class ParamsT> CppAD::vector<Type> loglikSeq(
const Rcpp::NumericVector Sstrata,
const CppAD::vector<Type>& gamma,
const PackedParams<ParamsT>& parameters, // gamma is ignored
const Data& data, 
const Config& cfg
){
  const size_t end = Sstrata.length();
  CppAD::vector<Type> result(1);
  result[0] = Type(0);

  for(size_t Diter=0; Diter < end; ++Diter ) {
    size_t Dstrata = Sstrata[Diter];
    auto etaHere = compute_eta_for_stratum<Type, ParamsT>(
      Dstrata, data, gamma, parameters.beta);

    auto contrib = accumulate_contrib_for_stratum<Type, ParamsT>(
      Dstrata, data, etaHere, parameters, cfg
      );
    result[0] += contrib[0];
  }
  return(result);
}




#endif
