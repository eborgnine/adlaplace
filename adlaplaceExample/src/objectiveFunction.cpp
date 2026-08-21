#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <cmath>
#include <utility>

#include "adlaplace/density_data.hpp"
#include "adlaplace/eta.hpp"
#include "adlaplace/rviews.hpp"
#include "adlaplace/atomics.hpp"

namespace {

std::pair<std::size_t, std::size_t> obs_range(
  const Config& config,
  std::size_t ny,
  std::size_t Dgroup) {

  const bool have_shards = config.obs_groups.ncol() > 0;
  std::size_t startP = 0;
  std::size_t endP = 0;
  if (have_shards) {
    startP = static_cast<std::size_t>(config.obs_groups.p[Dgroup]);
    endP = static_cast<std::size_t>(config.obs_groups.p[Dgroup + 1]);
  } else if (Dgroup == 0) {
    endP = ny;
  }
  return {startP, endP};
}

}  // namespace

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const density_data& model,
  const Config& config,
  const size_t Dgroup)
{

  CppAD::AD<double> result = 0.0;

  const std::size_t omega_global = model.theta_index(0);
  const std::size_t omega_row = model.theta_row(0);
  const std::size_t alpha_index = model.theta_index(1);
  CppAD::AD<double> omega_in = x[omega_global];
  CppAD::AD<double> omega =
    transform_theta_at(config, omega_row) ? CppAD::exp(omega_in) : omega_in;
  CppAD::AD<double> omega_sqrt2 = omega * SQRTTWO;
  CppAD::AD<double> alpha = x[alpha_index];

  const size_t ny = model.y.size();
  const auto range = obs_range(config, ny, Dgroup);
  const bool have_shards = config.obs_groups.ncol() > 0;

  for (size_t DI = range.first; DI < range.second; ++DI) {
    const size_t Dobs = have_shards ? static_cast<size_t>(config.obs_groups.i[DI]) : DI;

    const CppAD::AD<double> eta = eta_at(x, model, Dobs);
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
  const density_data& model,
  const Config& config)
{

  const std::size_t omega_global = model.theta_index(0);
  const std::size_t omega_row = model.theta_row(0);
  CppAD::AD<double> omega_in = x[omega_global];
  CppAD::AD<double> log_omega =
    transform_theta_at(config, omega_row) ? omega_in : CppAD::log(omega_in);

  const double ny = static_cast<double>(model.y.size());
  CppAD::AD<double> logDens = static_cast<double>(ny) * (CppAD::AD<double>(0.0) - log_omega)
    + static_cast<double>(ny) * (CppAD::AD<double>(0.0) - ONEHALFLOGTWOPI);

  CppAD::vector<CppAD::AD<double>> out_v(1);
  out_v[0] = logDens;
  return out_v;
}
