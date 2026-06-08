#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <cmath>

#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/math/constants.hpp"

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup)
{
  const Config config(config_list);

  CppAD::AD<double> result = 0.0;

  const std::size_t omega_index = model.theta_index(0);
  const std::size_t alpha_index = model.theta_index(1);
  CppAD::AD<double> omega_in = x[omega_index];
  CppAD::AD<double> omega = config.transform_theta ? CppAD::exp(omega_in) : omega_in;
  CppAD::AD<double> omega_sqrt2 = omega * SQRTTWO;
  CppAD::AD<double> alpha = x[alpha_index];

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

  if (config.verbose && Dgroup < 1) {
    Rcpp::Rcout << "logDensObs skewnormal_obs group " << Dgroup
      << " startP " << startP << " endP " << endP << " ny " << ny
      << " have_shards " << have_shards
      << " omega index " << omega_index << " alpha index " << alpha_index
      << " omega " << omega << " alpha " << alpha
      << " transform_theta " << config.transform_theta << "\n";
  }

  for (size_t DI = startP; DI < endP; ++DI) {
    const size_t Dobs = have_shards ? config.shards.i[DI] : DI;

    CppAD::AD<double> eta_fixed = 0.0;
    const size_t p0x = model.XTp.p[Dobs];
    const size_t p1x = model.XTp.p[Dobs + 1];
    for (size_t D = p0x; D < p1x; ++D) {
      eta_fixed += model.XTp.x[D] * x[model.XTp.i[D]];
    }

    CppAD::AD<double> eta_random = 0.0;
    const size_t p0a = model.ATp.p[Dobs];
    const size_t p1a = model.ATp.p[Dobs + 1];
    for (size_t D = p0a; D < p1a; ++D) {
      eta_random += model.ATp.x[D] * x[model.num_beta + model.ATp.i[D]];
    }

    CppAD::AD<double> eta = eta_fixed + eta_random;
    CppAD::AD<double> z = (CppAD::AD<double>(model.y[Dobs]) - eta) / omega_sqrt2;
    CppAD::AD<double> t = -alpha * z;

    result += -z * z + CppAD::log(CppAD::erfc(t)); 
  }

  CppAD::vector<CppAD::AD<double>> out(1);
  out[0] = result;
  return out;
}

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list)
{
  const Config config(config_list);

  const std::size_t omega_index = model.theta_index(0);
  CppAD::AD<double> omega_in = x[omega_index];
  CppAD::AD<double> log_omega = config.transform_theta ? omega_in : CppAD::log(omega_in);

  const double ny = static_cast<double>(model.y.size());
  CppAD::AD<double> logDens = static_cast<double>(ny) * (CppAD::AD<double>(0.0) - log_omega)
    + static_cast<double>(ny) * (CppAD::AD<double>(0.0) - ONEHALFLOGTWOPI);

  CppAD::vector<CppAD::AD<double>> out_v(1);
  out_v[0] = logDens;
  return out_v;
}
