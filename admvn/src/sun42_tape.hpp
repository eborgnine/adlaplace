#ifndef ADMVN_SUN42_TAPE_HPP
#define ADMVN_SUN42_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSun42D = 4;
constexpr std::size_t kSun42M = 2;
constexpr std::size_t kSun42NPar = 23;
constexpr std::size_t kSun42NDomain = kSun42NPar;
constexpr std::size_t kSun42JointD = kSun42D + kSun42M;
constexpr std::size_t kSun42HsNFree = kSun42JointD * (kSun42JointD - 1) / 2;

enum class Sun42ParMap {
  kBlockChol = 0,       // unimplemented for rectangular SUN(4,2)
  kHyperspherical = 1
};

struct Sun42Result {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct Sun42AdfunPack {
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

struct Sun42ObsShard {
  MvnTape p1_tape;
  std::vector<double> x;
  Sun42AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun42P2Shard {
  MvnTape p2_tape;
  Sun42AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun42TapeBundle {
  std::vector<Sun42ObsShard> shards;
  Sun42P2Shard p2;
  std::vector<double> weights;
  double weight_sum = 0.0;
  std::size_t n_obs = 0;
  int n_threads = 1;
  Sun42ParMap par_map = Sun42ParMap::kHyperspherical;
};

Sun42TapeBundle create_sun42_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1,
  const std::vector<double>& weights = {},
  Sun42ParMap par_map = Sun42ParMap::kHyperspherical);

Sun42Result eval_sun42_bundle(
  Sun42TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
