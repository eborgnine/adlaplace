#ifndef ADLAPLACE_BACKEND_HPP
#define ADLAPLACE_BACKEND_HPP

#include <cppad/cppad.hpp>
#include <Rcpp.h>
#include <Eigen/Sparse>
#include <vector>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/runtime/rviews.hpp"

// Stored copy of a Matrix dCHMsimpl object (symbolic shell from hessian_map()).
struct CholPattern {
	int n = 0;
	std::vector<int> p;
	std::vector<int> i;
	std::vector<int> nz;
	std::vector<int> nxt;
	std::vector<int> prv;
	std::vector<int> type;
	std::vector<int> colcount;
	std::vector<int> perm;
	// Inverse permutation for Eigen (pinv.indices()[orig] = j if perm[j] = orig); used with twistedBy(P).
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pinv;
	std::vector<int> dim;
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

// One AD shard (GroupPack) exposed through the C API.
using ad_vector = std::vector<adlaplace_adpack_handle*>;

using hessian_template = Eigen::SparseMatrix<int, Eigen::ColMajor, int>;


// Shard handles + Hessian templates/maps (filled from getAdFun() list via get_ad_groups()).
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

void ad_groups_attach_chol_pattern(ad_groups& groups, const Rcpp::List& ad_fun);

#endif
