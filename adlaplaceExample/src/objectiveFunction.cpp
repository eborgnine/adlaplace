#include "adlaplace/adlaplace.hpp"
#include "adlaplace/math/constants.hpp"
#include "adlaplace/logDens/random.hpp"

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const Data& data,
  const Config& config,
  const size_t Dgroup)
{
  CppAD::AD<double> result = 0.0;

  CppAD::AD<double> omega_in = x[config.theta_end - 2];
  CppAD::AD<double> omega = config.transform_theta ? CppAD::exp(omega_in) : omega_in;
  CppAD::AD<double> omega_sqrt2 = omega * CppAD::AD<double>(SQRTTWO);
  CppAD::AD<double> alpha = x[config.theta_end - 1];

  const bool have_groups = config.groups.ncol() > 0;
  const size_t startP = have_groups ? config.groups.p[Dgroup] : Dgroup;
  const size_t endP = have_groups ? config.groups.p[Dgroup + 1] : Dgroup + 1;

  for (size_t DI = startP; DI < endP; ++DI) {
    const size_t Dobs = have_groups ? config.groups.i[DI] : DI;

    CppAD::AD<double> eta_fixed = 0.0;
    for (size_t D = data.X.p[Dobs]; D < data.X.p[Dobs + 1]; ++D) {
      eta_fixed += data.X.x[D] * x[config.beta_begin + data.X.i[D]];
    }

    CppAD::AD<double> eta_random = 0.0;
    for (size_t D = data.A.p[Dobs]; D < data.A.p[Dobs + 1]; ++D) {
      eta_random += data.A.x[D] * x[config.gamma_begin + data.A.i[D]];
    }

    const CppAD::AD<double> eta = eta_fixed + eta_random;
    const CppAD::AD<double> z = (CppAD::AD<double>(data.y[Dobs]) - eta) / omega_sqrt2;
    const CppAD::AD<double> erfc_val = CppAD::erfc(-alpha * z);

    result += -z * z + CppAD::log(erfc_val);
  }

  CppAD::vector<CppAD::AD<double>> out(1);
  out[0] = result;
  return out;
}

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const Data& data,
  const Config& config)
{

  CppAD::AD<double> omega_in = x[config.theta_end - 2];
  CppAD::AD<double> log_omega = config.transform_theta ? omega_in : CppAD::log(omega_in);

  CppAD::vector<CppAD::AD<double>> out(1);
  out[0] = -CppAD::AD<double>(data.y.size()) * (log_omega + CppAD::AD<double>(ONEHALFLOGTWOPI));

  if(config.verbose) {
    Rcpp::Rcout << "logDensExtra " << out[0]  << " log omega " << log_omega << 
    " omega_in " << omega_in << " tr " << config.transform_theta << "\n";
  }


  return out;
}

#include "adlaplace/creators/adfun.hpp"
