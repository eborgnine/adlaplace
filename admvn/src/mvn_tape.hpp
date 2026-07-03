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
  unsigned int seed);

MvnResult eval_mvn_tape(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  bool inner = true);

}  // namespace admvn

#endif
