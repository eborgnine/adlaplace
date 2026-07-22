#ifndef ADMVN_MVN_ANALYTIC_HPP
#define ADMVN_MVN_ANALYTIC_HPP

#include "mvn_tape.hpp"

#include <vector>

namespace admvn {

// High-accuracy P(X <= h) for X ~ N(0, R) with R a correlation matrix,
// dimension 1..3 (univariate / bivariate / trivariate). Returns NaN if
// n is outside 1..3. Matches the sn/mnormt specialized low-d path in spirit:
// bivariate and trivariate use nested Gauss-Legendre + Phi / Phi_2.
double pmvn_cdf_std(
  const std::vector<double>& h,
  const std::vector<std::vector<double>>& R);

// P(X <= upper) for X ~ N(mean, Sigma) with all lower = -Inf and n <= 3.
// Returns NaN if the specialized path does not apply.
double pmvn_cdf_special(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  const std::vector<double>& lower);

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
