#ifndef ADMVN_MVN_ANALYTIC_HPP
#define ADMVN_MVN_ANALYTIC_HPP

#include "mvn_tape.hpp"

#include <vector>

namespace admvn {

// Analytic domain gradient of F = P(X <= upper), X ~ N(mean, Sigma),
// for lower = -Inf and n <= 3. Domain layout matches pack_domain with
// Sigma_ij = scale_i (L L^T)_ij scale_j (L = genz.ch, original index order).
// Returns empty if analytic path does not apply.
std::vector<double> analytic_mvn_domain_grad(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz,
  const std::vector<double>& lower);

// Analytic ∂F/∂upper and ∂F/∂mean from a covariance matrix in the same
// index order as upper/mean (n <= 3, lower = -Inf).
bool analytic_mvn_upper_mean_grad(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  const std::vector<double>& lower,
  std::vector<double>& grad_upper,
  std::vector<double>& grad_mean);

}  // namespace admvn

#endif
