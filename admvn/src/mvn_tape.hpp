#ifndef ADMVN_MVN_TAPE_HPP
#define ADMVN_MVN_TAPE_HPP

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

struct GenzPack {
  std::vector<double> scale;
  std::vector<std::vector<double>> ch;
};

struct MvnResult {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
  std::vector<double> gradient_mean;
};

struct MvnTape {
  std::size_t n = 0;
  std::size_t n_domain = 0;
  std::size_t n_vech = 0;
  std::vector<int> perm;
  std::vector<double> lower;
  std::vector<std::vector<std::vector<double>>> qmc_w;
  std::size_t n_points = 0;
  std::size_t n_shifts = 0;
  // When true, skip CppAD tape build; values use specialized d<=3 CDF
  // (fallback QMC) and gradients use analytic formulas (lower = -Inf).
  bool value_only = false;

  CppAD::ADFun<double> fun;
  CppAD::sparse_jac_work work_grad;
  CppAD::sparse_jac_work work_inner_grad;
  CppAD::sparse_hes_work work_hess;
  CppAD::sparse_hes_work work_inner_hess;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_grad;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_grad_inner;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_hessian;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_hessian_inner;
  CppAD::sparse_rc<CppAD::vector<size_t>> unused_pattern;
  CppAD::vector<double> w;
  CppAD::vector<double> x_seed;
};

std::size_t vech_size(std::size_t n);
std::size_t vech_index(std::size_t i, std::size_t j);
std::size_t domain_size(std::size_t n);

GenzPack pack_genz_ch(
  const std::vector<std::vector<double>>& sigma,
  const std::vector<int>& perm);

std::vector<double> pack_domain(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz);

MvnTape create_mvn_tape(
  const std::vector<double>& lower,
  const std::vector<double>& upper_seed,
  const std::vector<double>& mean_seed,
  const std::vector<std::vector<double>>& sigma_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  bool value_only = false);

// Compiled double QMC value (+ optional Monte Carlo error estimate).
double eval_mvn_value_double(
  const MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz,
  double* error_out);

MvnResult eval_mvn_tape(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  bool inner = true);

std::vector<double> eval_mvn_domain_grad(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma);

// Prefer analytic orthant gradients (n<=3, lower=-Inf); fall back to AD-QMC.
std::vector<double> eval_mvn_domain_grad_auto(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz);

}  // namespace admvn

#endif
