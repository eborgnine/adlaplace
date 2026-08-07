#ifndef ADMVN_SUN22_TAPE_HPP
#define ADMVN_SUN22_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSun22D = 2;
constexpr std::size_t kSun22M = 2;
constexpr std::size_t kSun22NPar = 10;
constexpr std::size_t kSun22NDomain = kSun22NPar;
constexpr std::size_t kSun22JointD = kSun22D + kSun22M;
constexpr std::size_t kSun22HsNFree = kSun22JointD * (kSun22JointD - 1) / 2;

enum class Sun22ParMap {
  kBlockChol = 0,
  kHyperspherical = 1
};

struct Sun22Result {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct Sun22AdfunPack {
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

struct Sun22ObsShard {
  MvnTape p1_tape;
  std::vector<double> x;
  Sun22AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun22P2Shard {
  MvnTape p2_tape;
  Sun22AdfunPack adfun;
  std::size_t pmvn_call_id = 0;
};

struct Sun22TapeBundle {
  std::vector<Sun22ObsShard> shards;
  Sun22P2Shard p2;
  std::vector<double> weights;
  double weight_sum = 0.0;
  std::size_t n_obs = 0;
  int n_threads = 1;
  Sun22ParMap par_map = Sun22ParMap::kBlockChol;
};

Sun22TapeBundle create_sun22_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1,
  const std::vector<double>& weights = {},
  Sun22ParMap par_map = Sun22ParMap::kBlockChol);

Sun22Result eval_sun22_bundle(
  Sun22TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
