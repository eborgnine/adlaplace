#include "lgamma.hpp"
#include <Rcpp.h>


// Define the global atomic function instance
atomic_lgamma_ad lgamma_ad_atomic("lgamma_ad_atomic");

// Constructor implementation
atomic_lgamma_ad::atomic_lgamma_ad(const std::string& name)
    : CppAD::atomic_four<double>(name) {}


bool atomic_lgamma_ad::for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y
  )  {
    type_y[0] = type_x[0];
    return true;
  }


bool atomic_lgamma_ad::forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
  ) {
//    size_t n_order = order_up + 1;
//    assert(tx.size() >= n_order);
//    assert(ty.size() >= n_order);

    double x0 = tx[0];
    double d1;

    // Zero order (function value)
    if (order_low <= 0) {
      ty[0] = std::lgamma(x0);
    }

    // First order (first derivative)
    if (order_up >= 1) {
      d1 = R::psigamma(x0, 0); // digamma
      double x1 = tx[1];
      ty[1] = d1 * x1;
    }

    // Second order (second derivative)
    if (order_up >= 2) {
      double d2 = R::psigamma(x0, 1); // trigamma
      double x1 = tx[1];
      double x2 = tx[2];
      ty[2] = d1 * x2 + d2 * x1 * x1;
      // Note: Removed 0.5 factor - CppAD handles Taylor coefficients
    }

    return true;
  }

 bool atomic_lgamma_ad::reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
  ) {

    std::fill(px.begin(), px.end(), 0.0);

    double x0 = tx[0];
    double d1, d2;

    // Order 0
    if (order_up >= 0) {
      d1 = R::psigamma(x0, 0);
      px[0] = py[0] * d1;
    }

    // Order 1
    if (order_up >= 1) {
      d2 = R::psigamma(x0, 1);
      double x1 = tx[1];
      px[0] += py[1] * d2 * x1;
      px[1] = py[1] * d1;
    }

    // Order 2
    if (order_up >= 2) {
      double x1 = tx[1];
      double x2 = tx[2];
      px[0] += py[2] * (d2 * x1 * x1 + d1 * x2);
      px[1] += py[2] * 2 * d2 * x1;
      px[2] = py[2] * d1;
    }

    return true;
  }



// Template function implementation
template<class Type>
Type lgamma_ad(Type x) {
    CppAD::vector<Type> xin(1), yout(1);
    xin[0] = x;
    lgamma_ad_atomic(xin, yout);
    return yout[0];
}

// Explicit instantiations
//template double lgamma_ad<double>(double);
template CppAD::AD<double> lgamma_ad<CppAD::AD<double>>(CppAD::AD<double>);
