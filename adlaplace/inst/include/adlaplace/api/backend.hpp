#ifndef ADLAPLACE_BACKEND_HPP
#define ADLAPLACE_BACKEND_HPP

#include <cppad/cppad.hpp>
#include <Rcpp.h>
#include <Eigen/Sparse>
#include <vector>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/creators/rviews.hpp"

// Symbolic LDL from hessian_map() chol_inner_list (CSC patterns, 0-based perm / perm_inv).
struct CholPattern {
	int n = 0;
	std::vector<int> L1_p;
	std::vector<int> L1_i;
	std::vector<int> Linv_p;
	std::vector<int> Linv_i;
	std::vector<int> perm;
	std::vector<int> perm_inv;
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

  // Shard index and global parameter layout (for trace_hinv_t / LinvPtColumns).
  std::size_t shard_index = 0;
  std::size_t n_beta = 0;
  std::size_t n_gamma = 0;
};

static inline GroupPack* pack_ctx(void* vctx) {
  return static_cast<GroupPack*>(vctx);
}

// One AD shard (GroupPack) exposed through the C API.
using ad_vector = std::vector<adlaplace_adpack_handle*>;

using hessian_template = Eigen::SparseMatrix<int, Eigen::ColMajor, int>;


// Shard handles + Hessian templates/maps (filled from get_ad_fun() list via get_ad_groups()).
struct ad_groups {
  ad_vector fun;
  hessian_template hessian_outer;
  hessian_template hessian_inner;
  hessian_map_view map_outer;
  hessian_map_view map_inner;
  IntVecView sizes;
  CholPattern chol_pattern;
  bool hessians_attached = false;
};

#endif
