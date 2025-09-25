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
  // tx[k] = x^(k), ty[k] = y^(k), n = m = 1
  const double x0 = tx[0];

  // f^(n)(x0): f1=digamma, f2=trigamma, f3=polygamma(2), f4=polygamma(3)
  const double f1 = (order_up >= 1 ? R::psigamma(x0, 0) : 0.0);
  const double f2 = (order_up >= 2 ? R::psigamma(x0, 1) : 0.0);
  const double f3 = (order_up >= 3 ? R::psigamma(x0, 2) : 0.0);
  const double f4 = (order_up >= 4 ? R::psigamma(x0, 3) : 0.0);

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
    ty[2] = f1 * x2 + 0.5 * f2 * (x1 * x1);                  // <-- 1/2
  }
  if (order_up >= 3) {
    const double x1 = tx[1];
    const double x2 = tx[2];
    const double x3 = tx[3];
    ty[3] = f1 * x3 + f2 * (x1 * x2) + (1.0/6.0) * f3 * (x1 * x1 * x1); // <-- 1, 1/6
  }
  if (order_up >= 4) {
    const double x1 = tx[1];
    const double x2 = tx[2];
    const double x3 = tx[3];
    const double x4 = tx[4];
    ty[4] = f1 * x4
          + f2 * ( x1 * x3 + 0.5 * (x2 * x2) )               // <-- + 1/2 x2^2
          + f3 * ( 0.5 * (x1 * x1 * x2) )                    // <-- 1/2
          + (1.0/24.0) * f4 * (x1 * x1 * x1 * x1);           // <-- 1/24
  }
  return true;
}

 bool atomic_lgamma_ad::reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& /*ty*/,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
) {
    std::fill(px.begin(), px.end(), 0.0);

    const double x0 = tx[0];

    // f^(n)(x0) for n=1..4 (deriv[0]=f', deriv[1]=f'', etc.)
    const double f1 = R::psigamma(x0, 0); // digamma
    const double f2 = R::psigamma(x0, 1); // trigamma
    const double f3 = R::psigamma(x0, 2);
    const double f4 = R::psigamma(x0, 3);

    // Convenience
    const double x1 = (order_up >= 1 ? tx[1] : 0.0);
    const double x2 = (order_up >= 2 ? tx[2] : 0.0);
    const double x3 = (order_up >= 3 ? tx[3] : 0.0);

    // --- order 0 ---
    // y^(0) depends only on x^(0): dy0/dx0 = f'(x0)
    if (order_up >= 0) {
        px[0] += py[0] * f1;
    }

    // --- order 1 ---
    // y^(1) = f'(x0) x1  =>  ∂y1/∂x0 = f'' x1,  ∂y1/∂x1 = f'
    if (order_up >= 1) {
        px[0] += py[1] * (f2 * x1);
        px[1] += py[1] * f1;
    }

    // --- order 2 ---
    // y^(2) = f' x2 + (1/2) f'' x1^2
    // ∂y2/∂x0 = f'' x2 + (1/2) f''' x1^2
    // ∂y2/∂x1 = f'' x1
    // ∂y2/∂x2 = f'
    if (order_up >= 2) {
        px[0] += py[2] * (f2 * x2 + 0.5 * f3 * x1 * x1);
        px[1] += py[2] * (f2 * x1);
        px[2] += py[2] * f1;
    }

    // --- order 3 ---
    // y^(3) = f' x3 + f'' x1 x2 + (1/6) f''' x1^3
    // ∂y3/∂x0 = f'' x3 + f''' x1 x2 + (1/6) f'''' x1^3
    // ∂y3/∂x1 = f'' x2 + (1/2) f''' x1^2
    // ∂y3/∂x2 = f'' x1
    // ∂y3/∂x3 = f'
    if (order_up >= 3) {
        px[0] += py[3] * (f2 * x3 + f3 * x1 * x2 + (1.0/6.0) * f4 * x1 * x1 * x1);
        px[1] += py[3] * (f2 * x2 + 0.5 * f3 * x1 * x1);
        px[2] += py[3] * (f2 * x1);
        px[3] += py[3] * f1;
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


