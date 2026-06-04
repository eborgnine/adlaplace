#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <cmath>

#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/math/constants.hpp"
#include "adlaplace/api/density_registry.hpp"

namespace {

constexpr double INV_SQRT_PI = 0.564189583547756286948079451560772585844050629328998810574693962167784;

class atomic_log_erfc_ad : public CppAD::atomic_four<double> {
public:
  explicit atomic_log_erfc_ad(const char* name)
    : CppAD::atomic_four<double>(name) {}

  bool for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y) override {
    (void)call_id;
    type_y.resize(1);
    type_y[0] = type_x[0];
    return true;
  }

  bool rev_depend(
    size_t call_id,
    const CppAD::vector<bool>& ident_zero_x,
    CppAD::vector<bool>& depend_x,
    const CppAD::vector<bool>& depend_y) override {
    (void)call_id;
    (void)ident_zero_x;
    depend_x.resize(1);
    depend_x[0] = depend_y[0];
    return true;
  }

  bool jac_sparsity(
    size_t call_id,
    bool dependency,
    const CppAD::vector<bool>& ident_zero_x,
    const CppAD::vector<bool>& select_x,
    const CppAD::vector<bool>& select_y,
    CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
    (void)call_id;
    (void)dependency;
    (void)ident_zero_x;
    if (select_x.size() != 1 || select_y.size() != 1) {
      return false;
    }
    if (!select_x[0] || !select_y[0]) {
      pattern_out.resize(1, 1, 0);
      return true;
    }
    pattern_out.resize(1, 1, 1);
    pattern_out.set(0, 0, 0);
    return true;
  }

  bool forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty) override {
    (void)call_id;
    (void)select_y;
    ty.resize(order_up + 1);
    const double x0 = tx[0];
    const double erfc_x = std::erfc(x0);
    if (order_low <= 0) {
      ty[0] = std::log(erfc_x);
    }
    if (order_up >= 1 && order_low <= 1) {
      const double x1 = tx[1];
      ty[1] = x1 * (-2.0 * INV_SQRT_PI * std::exp(-x0 * x0) / erfc_x);
    }
    return true;
  }

  bool reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py) override {
    (void)call_id;
    (void)select_x;
    (void)ty;
    px.resize(order_up + 1);
    for (size_t k = 0; k < px.size(); ++k) {
      px[k] = 0.0;
    }
    const double x0 = tx[0];
    const double erfc_x = std::erfc(x0);
    const double deriv0 = -2.0 * INV_SQRT_PI * std::exp(-x0 * x0) / erfc_x;
    px[0] += py[0] * deriv0;
    return true;
  }
};

inline CppAD::AD<double> log_erfc_ad(const CppAD::AD<double>& x) {
  thread_local static atomic_log_erfc_ad op("log_erfc_ad");
  CppAD::vector<CppAD::AD<double>> ax(1);
  CppAD::vector<CppAD::AD<double>> ay(1);
  ax[0] = x;
  op(ax, ay);
  return ay[0];
}

}  // namespace

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
  CppAD::AD<double> omega_sqrt2 = omega * CppAD::AD<double>(SQRTTWO);
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

    const CppAD::AD<double> eta = eta_fixed + eta_random;
    const CppAD::AD<double> z = (CppAD::AD<double>(model.y[Dobs]) - eta) / omega_sqrt2;
    const CppAD::AD<double> t = alpha * z;

    result += -z * z + log_erfc_ad(-t);
  }

  CppAD::vector<CppAD::AD<double>> out(1);
//  out[0] = result;
out[0] = omega_in;
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

  double ny = static_cast<double>(model.y.size());
  CppAD::AD<double> ny_ad = CppAD::AD<double>(ny), onehundred = 100;
  CppAD::AD<double> log_omega_ny = log_omega * onehundred;//ny_ad;

  CppAD::AD<double> out = log_omega_ny;// - ny * ONEHALFLOGTWOPI;

  CppAD::vector<CppAD::AD<double>> out_v(1);
  out_v[0] = out;
  return out_v;
}

// [[Rcpp::export]]
void register_example_densities() {
  register_log_dens_obs("skewnormal_obs", &logDensObs);
  register_log_dens_single_data("skewnormal_extra", &logDensExtra);
}
