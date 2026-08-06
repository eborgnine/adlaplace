#ifndef ADMVN_SUN44_TAPE_HPP
#define ADMVN_SUN44_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSun44D = 4;
constexpr std::size_t kSun44M = 4;
constexpr std::size_t kSun44NPar = 36;
constexpr std::size_t kSun44NDomain = kSun44NPar;

struct Sun44Result {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct Sun44AdfunPack {
  CppAD::ADFun<double> fun;
  CppAD::sparse_jac_work work_grad;
  CppAD::sparse_hes_work work_hess;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_grad;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_grad_inner;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_hessian;
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>> pattern_hessian_inner;
  CppAD::sparse_jac_work work_inner_grad;
  CppAD::sparse_hes_work work_inner_hess;
  CppAD::sparse_rc<CppAD::vector<size_t>> unused_pattern;
  CppAD::vector<double> w;
  CppAD::vector<double> x_seed;
};

struct Sun44ObsShard {
  MvnTape p1_tape;
  std::vector<double> x;
  Sun44AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun44P2Shard {
  MvnTape p2_tape;
  Sun44AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun44TapeBundle {
  std::vector<Sun44ObsShard> shards;
  Sun44P2Shard p2;
  std::vector<double> weights;
  double weight_sum = 0.0;
  std::size_t n_obs = 0;
  int n_threads = 1;
};

Sun44TapeBundle create_sun44_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1,
  const std::vector<double>& weights = {});

Sun44Result eval_sun44_bundle(
  Sun44TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
