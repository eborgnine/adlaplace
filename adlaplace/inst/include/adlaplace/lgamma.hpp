#pragma once

#ifndef LGAMMA_AD_HPP
#define LGAMMA_AD_HPP


#include <cppad/cppad.hpp>
#include <Rcpp.h>
#include <string>
#include <cmath>
#include <algorithm> 

class atomic_lgamma_ad : public CppAD::atomic_four<double> {
public:
  explicit atomic_lgamma_ad(const std::string& name)
  : CppAD::atomic_four<double>(name) {}

  // for_type
  bool for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y
    ) override {
    type_y.resize(1);
    type_y[0] = type_x[0];
    return true;
  }

bool rev_depend(
  size_t                                           call_id,
  const CppAD::vector<bool>&                      ident_zero_x,
  CppAD::vector<bool>&                            depend_x,
  const CppAD::vector<bool>&                      depend_y
) override
{
  // n = 1, m = 1
  // If y0 is relevant downstream, then x0 is relevant.
  depend_x.resize(1);
  depend_x[0] = depend_y[0];   // nothing fancy, y0 always depends on x0

  return true;
}

bool jac_sparsity(
  size_t                                           call_id,
  bool                                             dependency,
  const CppAD::vector<bool>&                      ident_zero_x,
  const CppAD::vector<bool>&                      select_x,
  const CppAD::vector<bool>&                      select_y,
  CppAD::sparse_rc< CppAD::vector<size_t> >&      pattern_out
) override
{
  const size_t n = select_x.size();  // domain size
  const size_t m = select_y.size();  // range size
  if (n != 1 || m != 1)
    return false;

  // If we don’t care about x0/y0, or x0 is identically zero,
  // treat Jacobian as structurally zero for this query.
  bool use_x0 = select_x[0] && select_y[0] && !ident_zero_x[0];
  if (!use_x0) {
    pattern_out.resize(m, n, 0);   // 1 x 1 zero matrix
    return true;
  }

  // y0 depends on x0 -> one nonzero at (0,0) in the Jacobian.
  pattern_out.resize(m, n, 1);
  pattern_out.set(0, 0, 0);
  return true;
}

  // forward
  bool forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
    ) override {
    ty.resize(order_up + 1);

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


  // reverse
      bool reverse(
        size_t call_id,
        const CppAD::vector<bool>& select_x,
        size_t order_up,
        const CppAD::vector<double>& tx,
        const CppAD::vector<double>& ty,
        CppAD::vector<double>& px,
        const CppAD::vector<double>& py
        ) override  {

        px.resize(order_up + 1);
        for(size_t i = 0; i < px.size(); ++i) {
          px[i] = 0.0;
        }

        const double x0 = tx[0];

  if (x0 <= 0.0){
    px[0] += py[0] * 1e10;
    return true;
  }

    // --- order 0 ---
    const double f1 = R::psigamma(x0, 0); // digamma
    const double f2 = order_up >= 1?R::psigamma(x0, 1):0.0;
    const double f3 = order_up >= 2?R::psigamma(x0, 2):0.0;
    const double f4 = order_up >= 3?R::psigamma(x0, 3):0.0;

    px[0] += py[0] * f1;

if (order_up >= 1) {
    const double x1 = tx[1];                     // requires Forward(1)
    px[0] += py[1] * (f2 * x1);
    px[1] += py[1] * f1;
  }

  if (order_up >= 2) {
    const double x1 = tx[1];                     // requires Forward(1)
    const double x2 = tx[2];                     // requires Forward(2)
    px[0] += py[2] * (f2 * x2 + 0.5 * f3 * x1 * x1);
    px[1] += py[2] * (f2 * x1);
    px[2] += py[2] * f1;
  }

  if (order_up >= 3) {
    const double x1 = tx[1];
    const double x2 = tx[2];
    const double x3 = tx[3];
    px[0] += py[3] * (f2 * x3 + f3 * x1 * x2 + (1.0/6.0) * f4 * x1 * x1 * x1);
    px[1] += py[3] * (f2 * x2 + 0.5 * f3 * x1 * x1);
    px[2] += py[3] * (f2 * x1);
    px[3] += py[3] * f1;
  }


    return true;
  }

};



inline CppAD::AD<double> lgamma_ad(const CppAD::AD<double>& x) {
    static atomic_lgamma_ad op("lgamma_ad");  // thread-safe since C++11
    CppAD::vector< CppAD::AD<double> > ax(1), ay(1);
    ax[0] = x;
  op(ax, ay);   // calls your atomic_four op
  return ay[0];
}


// declare primary template
template<class T>
inline T lgamma_any(T x);

// specializations (MUST be in a header)
template<>
inline double lgamma_any<double>(double x) {
  return std::lgamma(x);
}

template<>
inline CppAD::AD<double> lgamma_any< CppAD::AD<double> >(CppAD::AD<double> x) {
  // call your atomic wrapper (header-only, function-local static)
  extern CppAD::AD<double> lgamma_ad(const CppAD::AD<double>&);
  return lgamma_ad(x);
}

// declare primary template
template<class T>
inline T exp_any(T x);

// specializations (MUST be in a header)
template<>
inline double exp_any<double>(double x) {
  return std::exp(x);
}

template<>
inline CppAD::AD<double> exp_any< CppAD::AD<double> >(CppAD::AD<double> x) {
  return CppAD::exp(x);
}


#endif
