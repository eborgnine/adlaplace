#ifndef ADMVN_SUN43_TAPE_HPP
#define ADMVN_SUN43_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSun43D = 4;
constexpr std::size_t kSun43M = 3;
constexpr std::size_t kSun43NPar = 29;
constexpr std::size_t kSun43NDomain = kSun43NPar;
constexpr std::size_t kSun43JointD = kSun43D + kSun43M;
constexpr std::size_t kSun43HsNFree = kSun43JointD * (kSun43JointD - 1) / 2;

enum class Sun43ParMap {
  kBlockChol = 0,       // unimplemented for rectangular SUN(4,3)
  kHyperspherical = 1
};

struct Sun43Result {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct Sun43AdfunPack {
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

struct Sun43ObsShard {
  MvnTape p1_tape;
  std::vector<double> x;
  Sun43AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun43P2Shard {
  MvnTape p2_tape;
  Sun43AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun43TapeBundle {
  std::vector<Sun43ObsShard> shards;
  Sun43P2Shard p2;
  std::vector<double> weights;
  double weight_sum = 0.0;
  std::size_t n_obs = 0;
  int n_threads = 1;
  Sun43ParMap par_map = Sun43ParMap::kHyperspherical;
};

Sun43TapeBundle create_sun43_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1,
  const std::vector<double>& weights = {},
  Sun43ParMap par_map = Sun43ParMap::kHyperspherical);

Sun43Result eval_sun43_bundle(
  Sun43TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
