#ifndef LGAMMA_AD_HPP
#define LGAMMA_AD_HPP


#include <cppad/cppad.hpp>

class atomic_lgamma_ad : public CppAD::atomic_four<double> {
public:
  atomic_lgamma_ad(const std::string& name);

  // for_type
  bool for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y
  ) override;

  // forward
  bool forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
  ) override ;

  // reverse
  bool reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
  ) override;
};


// Declaration of the global atomic function instance
extern atomic_lgamma_ad lgamma_ad_atomic;

// Template function declaration
template<class Type>
Type lgamma_ad(Type x);


inline double lgamma_ad(double x) {
  return std::lgamma(x);
}

#endif
