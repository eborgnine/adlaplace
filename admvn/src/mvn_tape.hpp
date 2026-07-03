#ifndef ADMVN_MVN_TAPE_HPP
#define ADMVN_MVN_TAPE_HPP

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {

struct MvnSetup {
  std::size_t n = 0;
  std::vector<int> perm;
  std::vector<int> inv_perm;
  std::vector<double> mean;
  std::vector<double> scale;
  std::vector<double> as;
  std::vector<std::vector<double>> ch;
  std::vector<std::vector<std::vector<double>>> qmc_w;
  std::size_t n_points = 0;
  std::size_t n_shifts = 0;
};

struct MvnResult {
  double value = 0.0;
  double error = 0.0;
  std::vector<double> gradient;
  std::vector<std::vector<double>> hessian;
};

MvnSetup prepare_mvn_setup(
  const std::vector<double>& lower,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed);

MvnResult eval_mvn_tape(
  const MvnSetup& setup,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  CppAD::ADFun<double>* pre_taped = nullptr);

CppAD::ADFun<double> build_mvn_tape(const MvnSetup& setup);

}  // namespace admvn

#endif
