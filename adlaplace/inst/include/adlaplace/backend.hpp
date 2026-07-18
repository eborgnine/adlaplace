#ifndef ADLAPLACE_BACKEND_HPP
#define ADLAPLACE_BACKEND_HPP

#include <cppad/cppad.hpp>
#include <Rcpp.h>
#include <Eigen/SparseCore>
#include <cstddef>
#include <vector>

#include "adlaplace/ad_data.hpp"
#include "adlaplace/rviews.hpp"

// Bump when GroupPack / adlaplace_shard / ad_fun layout or the backend
// registration contract changes. packs_to_ad_fun stamps this into ad_fun;
// adlaplace checks it when consuming a backend-built handle.
#define ADLAPLACE_ABI_VERSION 1

// Symbolic LDL pattern from hessian_map()$chol_inner_list (L1, Linv, perm,
// half_H_inv, H_inv).
struct CholPattern {
  int n = 0;
  std::vector<int> L1_p;
  std::vector<int> L1_i;
  std::vector<int> Linv_p;
  std::vector<int> Linv_i;
  std::vector<int> half_H_inv_p;
  std::vector<int> half_H_inv_i;
  std::vector<int> H_inv_p;
  std::vector<int> H_inv_i;
  std::vector<int> trace_columns_p;
  std::vector<int> trace_columns_i;
  std::vector<int> perm;
  std::vector<int> perm_inv;
};

// Order-3 Forward/Reverse workspace for trace_hinv_t.
struct TraceWorkspace {
  CPPAD_TESTVECTOR(double) wthree;
  CPPAD_TESTVECTOR(double) direction_zeros;
  CPPAD_TESTVECTOR(double) direction;
};

struct GroupPack {
  CppAD::ADFun<double> fun;
  CppAD::sparse_jac_work work_grad;
  CppAD::sparse_hes_work work_hess;
  CppAD::sparse_jac_work work_inner_grad;
  CppAD::sparse_hes_work work_inner_hess;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad_inner;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian_inner;

  CPPAD_TESTVECTOR(double) w;
  TraceWorkspace trace;
  // Required but ignored output argument of sparse_jac_rev / sparse_hes.
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

  CPPAD_TESTVECTOR(double) x;

  std::size_t shard_index = 0;
  std::size_t owner_thread = 0;
  bool owner_thread_assigned = false;
  std::size_t n_beta = 0;
  std::size_t n_theta = 0;
};

inline GroupPack clone_group_pack(const GroupPack& src) {
  GroupPack dst;
  dst.fun = src.fun;
  dst.work_grad = src.work_grad;
  dst.work_hess = src.work_hess;
  dst.work_inner_grad = src.work_inner_grad;
  dst.work_inner_hess = src.work_inner_hess;
  dst.pattern_grad = src.pattern_grad;
  dst.pattern_grad_inner = src.pattern_grad_inner;
  dst.pattern_hessian = src.pattern_hessian;
  dst.pattern_hessian_inner = src.pattern_hessian_inner;
  dst.x = src.x;
  dst.w = src.w;
  dst.trace = src.trace;
  dst.unused_pattern = src.unused_pattern;
  dst.shard_index = src.shard_index;
  dst.owner_thread = src.owner_thread;
  dst.owner_thread_assigned = src.owner_thread_assigned;
  dst.n_beta = src.n_beta;
  dst.n_theta = src.n_theta;
  return dst;
}

struct adlaplace_shard;

using ShardFactory = adlaplace_shard* (*)(GroupPack&&);

struct adlaplace_shard {
  GroupPack pack;
  ShardFactory factory = nullptr;

  explicit adlaplace_shard(GroupPack&& p, ShardFactory f = nullptr)
    : pack(std::move(p)), factory(f) {}

  virtual ~adlaplace_shard() = default;

  virtual int f(const double* x, double* out_f) = 0;
  virtual int f_grad(const double* x, bool inner, double* out_f, double* out_grad) = 0;
  virtual int f_grad_hess(
    const double* x,
    bool inner,
    double* out_f,
    double* out_grad,
    double* out_hes,
    int* map) = 0;
  virtual int get_sizes(
    int* n_inner,
    int* n_outer,
    int* n_beta,
    int* n_theta,
    int* nnz_grad_inner,
    int* nnz_grad_outer,
    int* nnz_hes_inner,
    int* nnz_hes_outer) = 0;
  virtual int get_sparse_pattern(
    int* pattern_grad_inner,
    int* pattern_grad_outer,
    int* pattern_hes_inner_row,
    int* pattern_hes_inner_col,
    int* pattern_hes_outer_row,
    int* pattern_hes_outer_col) = 0;
  virtual int assign_memory() = 0;
  virtual int trace_hinv_t(
    const double* x,
    const int* LinvPt_p,
    const int* LinvPt_i,
    const double* LinvPt_x,
    std::size_t LinvPt_ncol,
    std::size_t LinvPt_p_len,
    std::size_t LinvPt_i_len,
    std::size_t LinvPt_x_len,
    const int* LinvPtColumns_p,
    const int* LinvPtColumns_i,
    std::size_t LinvPtColumns_p_len,
    std::size_t LinvPtColumns_i_len,
    double* out_trace) = 0;
  virtual adlaplace_shard* clone() const = 0;
};

using ad_vector = std::vector<adlaplace_shard*>;

using hessian_template = Eigen::SparseMatrix<int, Eigen::ColMajor, int>;

struct ad_fun {
  int abi_version = ADLAPLACE_ABI_VERSION;
  ad_vector fun;
  hessian_template hessian_outer;
  hessian_template hessian_inner;
  hessian_map_view map_outer;
  hessian_map_view map_inner;
  IntVecView sizes;
  CholPattern chol_pattern;
  bool hessians_attached = false;
  std::size_t configured_num_threads = 1;
  bool num_threads_configured = false;
};

using LogDensObsFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config,
  size_t Dgroup);

using LogDensSingleDataFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config);

using LogDensSingleRandomFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config);

#endif
