#ifndef ADMVN_SUN52_TAPE_HPP
#define ADMVN_SUN52_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSun52D = 5;
constexpr std::size_t kSun52M = 2;
constexpr std::size_t kSun52NPar = 31;
constexpr std::size_t kSun52NDomain = kSun52NPar;
constexpr std::size_t kSun52JointD = kSun52D + kSun52M;
constexpr std::size_t kSun52HsNFree = kSun52JointD * (kSun52JointD - 1) / 2;

enum class Sun52ParMap {
  kBlockChol = 0,       // unimplemented for rectangular SUN(5,2)
  kHyperspherical = 1
};

struct Sun52Result {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct Sun52AdfunPack {
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

struct Sun52ObsShard {
  MvnTape p1_tape;
  std::vector<double> x;
  Sun52AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun52P2Shard {
  MvnTape p2_tape;
  Sun52AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun52TapeBundle {
  std::vector<Sun52ObsShard> shards;
  Sun52P2Shard p2;
  std::vector<double> weights;
  double weight_sum = 0.0;
  std::size_t n_obs = 0;
  int n_threads = 1;
  Sun52ParMap par_map = Sun52ParMap::kHyperspherical;
};

Sun52TapeBundle create_sun52_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1,
  const std::vector<double>& weights = {},
  Sun52ParMap par_map = Sun52ParMap::kHyperspherical);

Sun52Result eval_sun52_bundle(
  Sun52TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
