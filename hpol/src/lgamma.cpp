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

    double deriv[4], tx1mult;

    // Zero order (function value)
    if (order_low <= 0) {
      ty[0] = std::lgamma(tx[0]);
    }

    if (order_up >= 1) { // 1st
      deriv[0] = R::psigamma(tx[0], 0); // digamma
      ty[1] = deriv[0] * tx[1];
    }

    if (order_up >= 2) { // 2nd
      deriv[1] = R::psigamma(tx[0], 1); // trigamma
      tx1mult = tx[1] * tx[1];
      ty[2] = deriv[0] * tx[2] + deriv[1] * tx1mult;
    }
    if (order_up >= 3) { // 3rd
      deriv[2] = R::psigamma(tx[0], 2); 
      tx1mult *= tx[1];
      ty[3] = deriv[0] * tx[3] + 3*deriv[1] * tx[1] * tx[2] +
        deriv[2] * tx1mult;
    }
  if (order_up >= 4) {
    deriv[3] = R::psigamma(tx[0], 3); // pentagamma
    tx1mult *= tx[1]; 
    ty[4] = deriv[0] * tx[4]
          + 4 * deriv[1] * tx[1] * tx[3]
          + 3 * deriv[1] * tx[2] * tx[2]
          + 6 * deriv[2] * tx[1] * tx[1] * tx[2]
          + deriv[3] * tx1mult;
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
    double deriv[4];
    double tx1sq, tx1cubed;


    // Order 0
    if (order_up >= 0) {
      deriv[0] = R::psigamma(x0, 0);
      px[0] = py[0] * deriv[0];
    }

    // Order 1
    if (order_up >= 1) {
      deriv[1] = R::psigamma(x0, 1);
      px[0] += py[1] * deriv[1] * tx[1];
      px[1] = py[1] * deriv[0];
    }

    // Order 2
    if (order_up >= 2) {
      deriv[2] = R::psigamma(x0, 2);      
      tx1sq = tx[1]*tx[1];
      px[0] += py[2] * (deriv[2] * tx1sq + deriv[1] * tx[2]);
      px[1] += py[2] * 2 * deriv[1] * tx[1];
      px[2] = py[2] * deriv[0];
    }

  // Order 3
    if (order_up >= 3) {
        tx1cubed = tx[1]*tx1sq;
        deriv[3] = R::psigamma(x0, 3);      
        px[0] += py[3] * 
          (deriv[1] * tx[3] + 
            3 * deriv[2] * tx[1] * tx[2] + 
            deriv[3] * tx1cubed);
        px[1] += py[3] * (
            3 * deriv[1] * tx[2] + 
            3 * deriv[2] * tx1sq);
        px[2] += py[3] * (3 * deriv[0] * tx[1]);
        px[3] += py[3] * deriv[0];
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
