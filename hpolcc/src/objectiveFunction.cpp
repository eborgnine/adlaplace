#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <cmath>

#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/creators/rviews.hpp"

static CppAD::AD<double> stable_logsumexp(
  const CppAD::vector<CppAD::AD<double>>& eta)
{
  CppAD::AD<double> max_eta = eta[0];
  for (size_t Deta = 1; Deta < eta.size(); ++Deta) {
    max_eta = CppAD::CondExpGt(eta[Deta], max_eta, eta[Deta], max_eta);
  }

  CppAD::AD<double> sumexp = 0.0;
  for (size_t Deta = 0; Deta < eta.size(); ++Deta) {
    sumexp += CppAD::exp(eta[Deta] - max_eta);
  }

  return CppAD::log(sumexp) + max_eta;
}

static CppAD::vector<CppAD::AD<double>> compute_eta_for_stratum(
  const size_t Dstrata,
  const ad_data& model,
  const CppAD::vector<CppAD::AD<double>>& x)
{
  const size_t startHere = model.elgm_matrix.p[Dstrata];
  const size_t endHere = model.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  CppAD::vector<CppAD::AD<double>> etaHere(NinStrata);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const size_t Deta = static_cast<size_t>(model.elgm_matrix.i[k]);

    CppAD::AD<double> accGamma = 0.0;
    CppAD::AD<double> accBeta = 0.0;

    const size_t p0x = model.XTp.p[Deta];
    const size_t p1x = model.XTp.p[Deta + 1];
    for (size_t t = p0x; t < p1x; ++t) {
      accBeta += model.XTp.x[t] * x[model.XTp.i[t]];
    }

    const size_t p0a = model.ATp.p[Deta];
    const size_t p1a = model.ATp.p[Deta + 1];
    for (size_t t = p0a; t < p1a; ++t) {
      accGamma += model.ATp.x[t] * x[model.num_beta + model.ATp.i[t]];
    }

    etaHere[j] = accBeta + accGamma;
  }

  return etaHere;
}

static CppAD::AD<double> accumulate_contrib_for_stratum(
  const size_t Dstrata,
  const ad_data& model,
  const CppAD::vector<CppAD::AD<double>>& etaHere,
  const CppAD::AD<double>& tauSq)
{
  const size_t startHere = model.elgm_matrix.p[Dstrata];
  const size_t endHere = model.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;

  const CppAD::AD<double> etaLogSum = stable_logsumexp(etaHere);
  CppAD::AD<double> contrib = 0.0;

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const size_t idx = static_cast<size_t>(model.elgm_matrix.i[k]);
    const double yhere = model.y[idx];

    if (yhere > 0) {
      const CppAD::AD<double> logMuKhere = etaHere[j] - etaLogSum;
      contrib += logMuKhere;

      if (yhere > 1) {
        const CppAD::AD<double> expMuKhere = CppAD::exp(logMuKhere);
        for (int i = 1; i < yhere; ++i) {
          const CppAD::AD<double> to_log = expMuKhere + static_cast<double>(i) * tauSq;
          contrib += CppAD::log(to_log);
        }
      }
    }
  }

  return contrib;
}

static inline CppAD::AD<double> tau_sq_from_x(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config)
{
  const std::size_t tau_index = model.theta_index(0);
  const CppAD::AD<double> tau_in = x[tau_index];
  const CppAD::AD<double> tau = config.transform_theta ? CppAD::exp(tau_in) : tau_in;
  return tau * tau;
}

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup)
{
  const Config config(config_list);
  CppAD::AD<double> result1 = 0.0;
  CppAD::vector<CppAD::AD<double>> result(1);

  const bool have_shards = config.shards.ncol() > 0;
  const size_t startP = have_shards ? config.shards.p[Dgroup] : Dgroup;
  const size_t endP = have_shards ? config.shards.p[Dgroup + 1] : Dgroup + 1;

  const CppAD::AD<double> tauSq = tau_sq_from_x(x, model, config);
  const bool verbose_here = config.verbose && (Dgroup < 1);

  if (verbose_here) {
    Rcpp::Rcout << "group " << Dgroup << " startP " << startP << " endP " << endP << "\n";
  }

  for (size_t DstrataI = startP; DstrataI < endP; ++DstrataI) {
    const size_t Dstrata = have_shards ? config.shards.i[DstrataI] : DstrataI;
    const auto etaHere = compute_eta_for_stratum(Dstrata, model, x);
    const auto contrib = accumulate_contrib_for_stratum(Dstrata, model, etaHere, tauSq);

    if (verbose_here) {
      Rcpp::Rcout << "strata " << Dstrata << " contrib " << contrib << "\n";
    }

    result1 += contrib;
  }

  if (verbose_here) {
    Rcpp::Rcout << " result " << result1 << "\n";
  }

  result[0] = result1;
  return result;
}

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list)
{
  const Config config(config_list);
  const CppAD::AD<double> tauSq = tau_sq_from_x(x, model, config);

  const size_t Nstrata = static_cast<size_t>(model.elgm_matrix.ncol());

  CppAD::AD<double> contribLog1pjTauSq = 0.0;
  double contribLogYkFact = 0.0;
  double contribLogNfact = 0.0;

  for (size_t Dstrata = 0; Dstrata < Nstrata; ++Dstrata) {
    int sumYhere = 0;
    const size_t endHere = model.elgm_matrix.p[Dstrata + 1];

    for (size_t DobsI = model.elgm_matrix.p[Dstrata]; DobsI < endHere; ++DobsI) {
      const size_t Dobs = static_cast<size_t>(model.elgm_matrix.i[DobsI]);
      sumYhere += static_cast<int>(model.y[Dobs]);
      contribLogYkFact += std::lgamma(1.0 + model.y[Dobs]);
    }

    contribLogNfact += std::lgamma(1.0 + sumYhere);
    const size_t sumYhereS = static_cast<size_t>(sumYhere);
    for (size_t Dsum = 1; Dsum < sumYhereS; ++Dsum) {
      contribLog1pjTauSq += CppAD::log1p(static_cast<double>(Dsum) * tauSq);
    }
  }

  const CppAD::AD<double> contrib = contribLogNfact - contribLogYkFact - contribLog1pjTauSq;

  if (config.verbose) {
    Rcpp::Rcout << "logDensExtra " << contrib << " tauSq " << tauSq
      << " contribLogNfact " << contribLogNfact
      << " contribLogYkFact " << contribLogYkFact
      << " contribLog1pjTauSq " << contribLog1pjTauSq << "\n";
  }

  CppAD::vector<CppAD::AD<double>> out(1);
  out[0] = contrib;
  return out;
}
