#ifndef ADLAPLACE_BACKEND_HPP
#define ADLAPLACE_BACKEND_HPP

#include <cppad/cppad.hpp>
#include <vector>

// Aggregated Hessian template + shard maps (built in R via hessianMap).
// Not stored on the AD handle; kept here for Rcpp conversion helpers.
struct HessianPack {
  std::vector<int> hessian_p;
  std::vector<int> hessian_i;
  std::vector<int> dim;

  std::vector<int> map_p;
  std::vector<int> map_local;
  std::vector<int> map_global;
};

struct GroupPack {
  CppAD::ADFun<double>              fun;
  CppAD::sparse_jac_work            work_grad;
  CppAD::sparse_hes_work            work_hess;
  CppAD::sparse_jac_work            work_inner_grad;
  CppAD::sparse_hes_work            work_inner_hess;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad_inner;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian_inner;

  CPPAD_TESTVECTOR(double) w, wthree, direction_zeros, direction;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

  CPPAD_TESTVECTOR(double) x;
};

using AdGroups = std::vector<GroupPack>;

#endif
