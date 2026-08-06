#pragma once

#include <cppad/cppad.hpp>

namespace admvn {

/// Floor an orthant probability before taking log so reverse-mode AD is zero
/// when p is at/below the floor (avoids 1/(p+eps) blow-ups when p underflows).
inline CppAD::AD<double> log_prob_floored(
  const CppAD::AD<double>& p,
  double eps = 1e-16) {
  using AD = CppAD::AD<double>;
  const AD e = AD(eps);
  const AD p_safe = CppAD::CondExpGt(p, e, p, e);
  return CppAD::log(p_safe);
}

} // namespace admvn
