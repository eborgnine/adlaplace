#ifndef ADLAPLACE_LOG_TOOLS_HPP
#define ADLAPLACE_LOG_TOOLS_HPP

// Numerically stable log helpers for AD densities: softplus, log(e^a + e^b),
// log-sum-exp, and log(e^a - e^b). Header-only (inline), so single-file
// sourceCpp models and backend packages can use them without linking
// adlaplace.so. Include directly; sourcecpp.hpp and atomics.hpp do not
// include this header.

#include <cppad/cppad.hpp>

namespace adlaplace {

// log(1 + exp(z)), stable for large |z|: exp(z) would overflow for z > 700.
inline CppAD::AD<double> softplus(const CppAD::AD<double>& z) {
  return CppAD::CondExpGt(
    z, CppAD::AD<double>(0.0),
    z + CppAD::log1p(CppAD::exp(-z)),
    CppAD::log1p(CppAD::exp(z)));
}

// log(exp(a) + exp(b)), stable for large negative inputs. CondExpGt is
// nondifferentiable on a == b, but the log-sum-exp correction restores the
// softmax derivative there (including the tie).
inline CppAD::AD<double> logaddexp(
  const CppAD::AD<double>& a,
  const CppAD::AD<double>& b) {

  const CppAD::AD<double> m = CppAD::CondExpGt(a, b, a, b);
  return m + CppAD::log(CppAD::exp(a - m) + CppAD::exp(b - m));
}

// log(sum(exp(eta))). eta must be non-empty.
inline CppAD::AD<double> logsumexp(
  const CppAD::vector<CppAD::AD<double>>& eta) {

  CppAD::AD<double> max_eta = eta[0];
  for (std::size_t Deta = 1; Deta < eta.size(); ++Deta) {
    max_eta = CppAD::CondExpGt(eta[Deta], max_eta, eta[Deta], max_eta);
  }

  CppAD::AD<double> sumexp = 0.0;
  for (std::size_t Deta = 0; Deta < eta.size(); ++Deta) {
    sumexp += CppAD::exp(eta[Deta] - max_eta);
  }

  return CppAD::log(sumexp) + max_eta;
}

// log(exp(a) - exp(b)). Requires a > b; a == b gives -Inf and a < b gives NaN.
inline CppAD::AD<double> logsubexp(
  const CppAD::AD<double>& a,
  const CppAD::AD<double>& b) {

  return a + CppAD::log1p(-CppAD::exp(b - a));
}

}  // namespace adlaplace

#endif
