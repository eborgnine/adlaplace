#ifndef ADMVN_GENZ_QMC_HPP
#define ADMVN_GENZ_QMC_HPP

#include "qnorm_atomic.hpp"

#include <cmath>
#include <vector>

namespace admvn {
namespace genz {

constexpr double kNine = 9.0;

template <typename ADType>
using ChMatrix = std::vector<std::vector<ADType>>;

template <typename ADType>
ADType qmc_integrand(
  const std::vector<ADType>& bs,
  const std::vector<ADType>& as,
  const ChMatrix<ADType>& ch,
  const std::vector<double>& w) {

  const std::size_t n = bs.size();
  const ADType ct0 = ch[0][0];
  const ADType ai0 = as[0];
  const ADType bi0 = bs[0];

  ADType c1 = ADType(0.0);
  if (ai0 > ADType(-kNine) * ct0) {
    if (ai0 < ADType(kNine) * ct0) {
      c1 = detail::pnorm_ad(ai0 / ct0);
    } else {
      c1 = ADType(1.0);
    }
  }

  ADType d1 = ADType(0.0);
  if (bi0 > ADType(-kNine) * ct0) {
    if (bi0 < ADType(kNine) * ct0) {
      d1 = detail::pnorm_ad(bi0 / ct0);
    } else {
      d1 = ADType(1.0);
    }
  }

  ADType c_val = c1;
  ADType dc_val = d1 - c1;
  ADType pv = dc_val;
  std::vector<ADType> y(n - 1, ADType(0.0));

  for (std::size_t i = 1; i < n; ++i) {
    const ADType x = ADType(w[i - 1]);
    y[i - 1] = detail::qnorm_ad(c_val + x * dc_val);

    ADType s = ADType(0.0);
    for (std::size_t k = 0; k < i; ++k) {
      s += ch[i][k] * y[k];
    }

    const ADType ct = ch[i][i];
    const ADType aicnt = as[i] - s;
    const ADType bicnt = bs[i] - s;

    ADType ci = ADType(1.0);
    ADType di = ADType(1.0);

    if (aicnt < ADType(-kNine) * ct) {
      ci = ADType(0.0);
    } else if (CppAD::abs(aicnt) < ADType(kNine) * ct) {
      ci = detail::pnorm_ad(aicnt / ct);
    }

    if (bicnt < ADType(-kNine) * ct) {
      di = ADType(0.0);
    } else if (CppAD::abs(bicnt) < ADType(kNine) * ct) {
      di = detail::pnorm_ad(bicnt / ct);
    }

    dc_val = di - ci;
    pv *= dc_val;
    c_val = ci;
  }

  return pv;
}

}  // namespace genz
}  // namespace admvn

#endif
