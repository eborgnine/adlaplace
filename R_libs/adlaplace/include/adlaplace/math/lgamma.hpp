#ifndef ADLAPLACE_MATH_LGAMMA_HPP
#define ADLAPLACE_MATH_LGAMMA_HPP

#include <cppad/cppad.hpp>
#include <cmath>
#include <string>

#ifdef USEBOOST
#include <boost/math/special_functions/polygamma.hpp>
#else
#include <Rcpp.h>
#include <Rmath.h>
#endif

#if defined(__GNUC__)
#define ADLAPLACE_MATH_EXPORT __attribute__((visibility("default")))
#else
#define ADLAPLACE_MATH_EXPORT
#endif

namespace adlaplace {
namespace detail {

constexpr double ONESIXTH = 1.0 / 6.0;
constexpr double ONEOVERTWENTYFOUR = 1.0 / 24.0;

inline double polygamma_threadsafe(const double x, const unsigned deriv) {
#ifdef USEBOOST
  return boost::math::polygamma(deriv, x);
#else
  return R::psigamma(x, deriv);
#endif
}

class atomic_lgamma_ad : public CppAD::atomic_four<double> {
public:
  explicit atomic_lgamma_ad(const std::string& name)
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
    const size_t n = select_x.size();
    const size_t m = select_y.size();
    if (n != 1 || m != 1) {
      return false;
    }
    const bool use_x0 = select_x[0] && select_y[0] && !ident_zero_x[0];
    if (!use_x0) {
      pattern_out.resize(m, n, 0);
      return true;
    }
    pattern_out.resize(m, n, 1);
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
    const double f1 = (order_up >= 1 ? polygamma_threadsafe(x0, 0U) : 0.0);
    const double f2 = (order_up >= 2 ? polygamma_threadsafe(x0, 1U) : 0.0);
    const double f3 = (order_up >= 3 ? polygamma_threadsafe(x0, 2U) : 0.0);
    const double f4 = (order_up >= 4 ? polygamma_threadsafe(x0, 3U) : 0.0);

    if (order_low <= 0) {
      ty[0] = std::lgamma(x0);
    }
    if (order_up >= 1 && order_low <= 1) {
      const double x1 = tx[1];
      ty[1] = f1 * x1;
    }
    if (order_up >= 2) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      ty[2] = f1 * x2 + 0.5 * f2 * (x1 * x1);
    }
    if (order_up >= 3) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      const double x3 = tx[3];
      ty[3] = f1 * x3 + f2 * (x1 * x2) + ONESIXTH * f3 * (x1 * x1 * x1);
    }
    if (order_up >= 4) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      const double x3 = tx[3];
      const double x4 = tx[4];
      ty[4] = f1 * x4
        + f2 * (x1 * x3 + 0.5 * (x2 * x2))
        + f3 * (0.5 * (x1 * x1 * x2))
        + ONEOVERTWENTYFOUR * f4 * (x1 * x1 * x1 * x1);
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
    for (size_t i = 0; i < px.size(); ++i) {
      px[i] = 0.0;
    }

    const double x0 = tx[0];
    if (x0 <= 0.0) {
      px[0] += py[0] * 1e10;
      return true;
    }

    const double f1 = polygamma_threadsafe(x0, 0U);
    const double f2 = order_up >= 1 ? polygamma_threadsafe(x0, 1U) : 0.0;
    const double f3 = order_up >= 2 ? polygamma_threadsafe(x0, 2U) : 0.0;
    const double f4 = order_up >= 3 ? polygamma_threadsafe(x0, 3U) : 0.0;

    px[0] += py[0] * f1;

    if (order_up >= 1) {
      const double x1 = tx[1];
      px[0] += py[1] * (f2 * x1);
      px[1] += py[1] * f1;
    }

    if (order_up >= 2) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      px[0] += py[2] * (f2 * x2 + 0.5 * f3 * x1 * x1);
      px[1] += py[2] * (f2 * x1);
      px[2] += py[2] * f1;
    }

    if (order_up >= 3) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      const double x3 = tx[3];
      px[0] += py[3] * (f2 * x3 + f3 * x1 * x2 + ONESIXTH * f4 * x1 * x1 * x1);
      px[1] += py[3] * (f2 * x2 + 0.5 * f3 * x1 * x1);
      px[2] += py[3] * (f2 * x1);
      px[3] += py[3] * f1;
    }

    return true;
  }

  bool hes_sparsity(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    const CppAD::vector<bool>& select_y,
    CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
    (void)call_id;
    const size_t n = select_x.size();
    const size_t m = select_y.size();
    if (n != 1 || m != 1) {
      return false;
    }
    const bool use_x0 = select_x[0] && select_y[0];
    if (!use_x0) {
      pattern_out.resize(n, n, 0);
      return true;
    }
    pattern_out.resize(n, n, 1);
    pattern_out.set(0, 0, 0);
    return true;
  }
};

}  // namespace detail
}  // namespace adlaplace

#ifdef ADLAPLACE_MATH_LGAMMA_DEFINE

namespace adlaplace {
namespace detail {

inline atomic_lgamma_ad& lgamma_atomic_instance() {
  static atomic_lgamma_ad op("lgamma_ad");
  return op;
}

}  // namespace detail
}  // namespace adlaplace

ADLAPLACE_MATH_EXPORT void adlaplace_init_lgamma_atomic() {
  (void)adlaplace::detail::lgamma_atomic_instance();
}

ADLAPLACE_MATH_EXPORT CppAD::AD<double> lgamma_ad(const CppAD::AD<double>& x) {
  CppAD::vector<CppAD::AD<double>> ax(1);
  CppAD::vector<CppAD::AD<double>> ay(1);
  ax[0] = x;
  adlaplace::detail::lgamma_atomic_instance()(ax, ay);
  return ay[0];
}

#else

ADLAPLACE_MATH_EXPORT void adlaplace_init_lgamma_atomic();
ADLAPLACE_MATH_EXPORT CppAD::AD<double> lgamma_ad(const CppAD::AD<double>& x);

#endif

#endif
