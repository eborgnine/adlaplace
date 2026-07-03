#ifndef ADMVN_QNORM_ATOMIC_HPP
#define ADMVN_QNORM_ATOMIC_HPP

#include <cppad/cppad.hpp>
#include <cmath>
#include <string>

#include <Rcpp.h>
#include <Rmath.h>

namespace admvn {
namespace detail {

inline double dnorm_std(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

class atomic_qnorm : public CppAD::atomic_four<double> {
public:
  explicit atomic_qnorm(const std::string& name)
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

  bool hes_sparsity(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    const CppAD::vector<bool>& select_y,
    CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
    (void)call_id;
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

    const double p0 = tx[0];
    const double q0 = R::qnorm(p0, 0.0, 1.0, 1, 0);
    const double phi0 = dnorm_std(q0);
    const double qp0 = 1.0 / phi0;
    const double qpp0 = q0 / (phi0 * phi0);

    if (order_low <= 0) {
      ty[0] = q0;
    }
    if (order_up >= 1 && order_low <= 1) {
      ty[1] = qp0 * tx[1];
    }
    if (order_up >= 2) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      ty[2] = qpp0 * (x1 * x1) + qp0 * x2;
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

    const double p0 = tx[0];
    const double q0 = R::qnorm(p0, 0.0, 1.0, 1, 0);
    const double phi0 = dnorm_std(q0);
    const double qp0 = 1.0 / phi0;
    const double qpp0 = q0 / (phi0 * phi0);

    px[0] += py[0] * qp0;

    if (order_up >= 1) {
      const double x1 = tx[1];
      px[0] += py[1] * (qpp0 * x1);
      px[1] += py[1] * qp0;
    }

    if (order_up >= 2) {
      const double x1 = tx[1];
      const double x2 = tx[2];
      px[0] += py[2] * (qpp0 * x2 + qpp0 * x1);  // d/dp of qp0*x1^2 + qp0*x2
      px[1] += py[2] * (2.0 * qpp0 * x1);
      px[2] += py[2] * qp0;
    }

    return true;
  }
};

inline atomic_qnorm& qnorm_atomic_instance() {
  static atomic_qnorm op("admvn_qnorm");
  return op;
}

inline void init_qnorm_atomic() {
  (void)qnorm_atomic_instance();
}

template <typename ADType>
inline ADType qnorm_ad(const ADType& p_raw) {
  const ADType eps = ADType(1e-12);
  const ADType one = ADType(1.0);
  const ADType p = CppAD::CondExpGt(
    p_raw, one - eps, one - eps,
    CppAD::CondExpLt(p_raw, eps, eps, p_raw));
  CppAD::vector<ADType> ax(1);
  CppAD::vector<ADType> ay(1);
  ax[0] = p;
  qnorm_atomic_instance()(ax, ay);
  return ay[0];
}

template <typename ADType>
inline ADType pnorm_ad(const ADType& x) {
  static const double inv_sqrt2 = 0.70710678118654752440;
  return ADType(0.5) * CppAD::erfc(-x * inv_sqrt2);
}

}  // namespace detail
}  // namespace admvn

#endif
