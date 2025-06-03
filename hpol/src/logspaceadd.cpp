#include "logspaceadd.hpp"
#include <Rcpp.h>

// Define the global atomic function instance
atomic_logspace_add logspace_add_atomic("logspace_add_atomic");

// Constructor implementation
atomic_logspace_add::atomic_logspace_add(const std::string& name)
    : CppAD::atomic_four<double>(name) {}


bool atomic_logspace_add::for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y
  )  {
    assert(type_x.size() == 2);
    assert(type_y.size() == 1);
    type_y[0] = std::max(type_x[0], type_x[1]);
    return true;
  }

bool atomic_logspace_add::forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
  )  {
    size_t n_order = order_up + 1;
    assert(tx.size() >= 2 * n_order);
    assert(ty.size() >= 1 * n_order);

    // Zeroth-order values
    double x0 = tx[0];
    double x1 = tx[1];
 //   double x2 = tx[2];


    double y0 = tx[n_order];
    double y1 = tx[n_order + 1];
//    double y2 = tx[n_order + 2];

    int xIsSmaller = x0 < y0;
    double diff0 = (xIsSmaller ? 
      y0-x0 :
      x0-y0 
      );
    double log1pexpdiff0 = std::log1p(std::exp(-diff0));
    double res, logDx, logDy;
     if(xIsSmaller) {
      res = y0 + log1pexpdiff0; // log(exp(x) + exp(y))
      logDx = - diff0 - log1pexpdiff0;
      logDy = - log1pexpdiff0;
    } else {
      res = x0 + log1pexpdiff0;
      logDx = - log1pexpdiff0;
      logDy = - diff0 - log1pexpdiff0;
    }
    double logD2 = x0+y0-2*res;
    double dx = std::exp(logDx), dy=std::exp(logDy), d2 = std::exp(logD2);



    // Zero order (function value)
    if (order_low <= 0) {
      ty[0] = res;
    }

    // First order (first derivatives)
    if (order_up >= 1) {
      ty[1] = dx * x1 + dy * y1;
    }

    // Second order (second derivatives)
    if (order_up >= 2) {
    double x1 = tx[2];  // x, order 1
    double y1 = tx[n_order + 2];  // y, order 1
    double x2 = tx[4];  // x, order 2
    double y2 = tx[n_order + 4];  // y, order 2

    // Function values at zero order (already computed: x0, y0, res, dx, dy, d2)
    // d2 = exp(x0 + y0 - 2 * res)

    // Second order Taylor coefficient for the output
    ty[2] = dx * x2 + dy * y2 + d2 * (x1 - y1) * (x1 - y1);
    }

    return true;
  }

bool atomic_logspace_add::reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
)  {
    size_t q = order_up + 1;
    assert(tx.size() >= 2 * q);
    assert(ty.size() >= 1 * q);
    assert(px.size() >= 2 * q);
    assert(py.size() >= 1 * q);

    std::fill(px.begin(), px.end(), 0.0);

    double x0 = tx[0];
    double x1 = tx[1];
//    double x2 = tx[2];

    double y0 = tx[q];
    double y1 = tx[q + 1];
//    double y2 = tx[q + 2];

    int xIsSmaller = x0 < y0;
    double diff = x1 - y1;
    double diff0 = (xIsSmaller ? 
      y0-x0 :
      x0-y0 
      );

    double log1pexpdiff0 = std::log1p(std::exp(-diff0));
    double res, logDx, logDy;
    if(xIsSmaller) {
      res = y0 + log1pexpdiff0; // log(exp(x) + exp(y))
      logDx = - diff0 - log1pexpdiff0;
      logDy = - log1pexpdiff0;
    } else {
      res = x0 + log1pexpdiff0;
      logDx = - log1pexpdiff0;
      logDy = - diff0 - log1pexpdiff0;
    }
    double logD2 = x0+y0-2*res;
    double dx = std::exp(logDx), dy=std::exp(logDy);
    double d2 = std::exp(logD2);


 
    // Order 0
    if (order_up >= 0) {
        px[0] = py[0] * dx;
        px[q] = py[0] * dy;
    }

    // Order 1
    if (order_up >= 1) {

        px[0] += py[1] * d2 * diff;   
        px[1] = py[1] * dx;

        px[q] += py[1] * d2 * (-diff);   
        px[q + 1] = py[1] * dy;
    }

    // Order 2
    if (order_up >= 2) {
        double x1 = tx[1];
        double y1 = tx[q + 1];
//        double x2 = tx[2];
//        double y2 = tx[q + 2];

        double dx1 = x1 - y1;
        // Second derivatives of logsumexp are:
        // d2x = d2 * (1 - 2*dx)
        // d2y = d2 * (1 - 2*dy)
        // d2xy = -d2

        // Output is z2 = dx * x2 + dy * y2 + d2 * (x1 - y1)^2

        // dL/dx2 = py[2] * dx
        px[2] += py[2] * dx;
        // dL/dy2 = py[2] * dy
        px[q + 2] += py[2] * dy;

        // dL/dx1 = py[2] * (2*d2*(x1 - y1))
        px[1] += py[2] * 2.0 * d2 * dx1;
        px[q + 1] += py[2] * 2.0 * d2 * (-dx1);

    }
    return true;
}

template<class Type>
Type logspace_add_ad(Type x, Type y) {
  CppAD::vector<Type> tx(2), ty(1);
  tx[0] = x;
  tx[1] = y;
  logspace_add_atomic(tx, ty);
  return ty[0];
}


//template double logspace_add_ad<double>(double, double);
template CppAD::AD<double> logspace_add_ad<CppAD::AD<double>>(CppAD::AD<double>, CppAD::AD<double>);







