#ifndef ADMVN_SUN_TAPE_HPP
#define ADMVN_SUN_TAPE_HPP

#include "mvn_tape.hpp"

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

constexpr std::size_t kSunD = 3;
constexpr std::size_t kSunM = 3;
constexpr std::size_t kSunNPar = 21;
constexpr std::size_t kSunNDomain = kSunNPar;

struct SunResult {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

struct SunAdfunPack {
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

struct SunObsShard {
  MvnTape p1_tape;
  std::array<double, kSunD> x{};
  SunAdfunPack adfun;
};

struct SunP2Shard {
  MvnTape p2_tape;
  SunAdfunPack adfun;
};

struct SunTapeBundle {
  std::vector<SunObsShard> shards;
  SunP2Shard p2;
  std::size_t n_obs = 0;
  int n_threads = 1;
};

SunTapeBundle create_sun_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads = 1);

SunResult eval_sun_bundle(
  SunTapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0,
  bool manage_parallel_scope = true);

}  // namespace admvn

#endif
