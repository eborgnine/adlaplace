#ifndef ADMVN_SUN32_TAPE_HPP
#define ADMVN_SUN32_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSun32D = 3;
constexpr std::size_t kSun32M = 2;
constexpr std::size_t kSun32NPar = 16;
constexpr std::size_t kSun32NDomain = kSun32NPar;
constexpr std::size_t kSun32JointD = kSun32D + kSun32M;
constexpr std::size_t kSun32HsNFree = kSun32JointD * (kSun32JointD - 1) / 2;

enum class Sun32ParMap {
  kBlockChol = 0,       // unimplemented for rectangular SUN(3,2)
  kHyperspherical = 1
};

struct Sun32Result {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct Sun32AdfunPack {
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

struct Sun32ObsShard {
  MvnTape p1_tape;
  std::vector<double> x;
  Sun32AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun32P2Shard {
  MvnTape p2_tape;
  Sun32AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun32TapeBundle {
  std::vector<Sun32ObsShard> shards;
  Sun32P2Shard p2;
  std::vector<double> weights;
  double weight_sum = 0.0;
  std::size_t n_obs = 0;
  int n_threads = 1;
  Sun32ParMap par_map = Sun32ParMap::kHyperspherical;
};

Sun32TapeBundle create_sun32_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1,
  const std::vector<double>& weights = {},
  Sun32ParMap par_map = Sun32ParMap::kHyperspherical);

Sun32Result eval_sun32_bundle(
  Sun32TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
