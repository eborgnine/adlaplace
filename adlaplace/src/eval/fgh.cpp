#include "adlaplace/eval/fgh.hpp"

#include <cppad/cppad.hpp>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/runtime/sizes.hpp"

int get_hessian(
  void* vctx,
  const bool* inner,
  const int** p,
  size_t* p_len,
  const int** i,
  size_t* i_len,
  int** map,
  size_t* map_len) {
  (void)vctx;
  (void)inner;
  (void)p;
  (void)p_len;
  (void)i;
  (void)i_len;
  (void)map;
  (void)map_len;
  return 1;
}

int eval_f(void* vctx, const double* x, double* out_f) {
  GroupPack& gp = *pack_ctx(vctx);
  const size_t Nparams = gp.x.size();
  const size_t Ndomain = gp.fun.Domain();
  const size_t Nrange = gp.fun.Range();
  if (Ndomain != Nparams) return 4;
  if (Nrange < 1) return 5;

  for (size_t D = 0; D < Nparams; D++) {
    gp.x[D] = x[D];
  }
  CppAD::vector<double> y = gp.fun.Forward(0, gp.x);
  if (y.size() < 1) return 6;
  *out_f += y[0];
  return 0;
}

int eval_grad(void* vctx, const double* x, const bool* inner, double* out_f, double* out_grad) {
  const bool innerv = *inner;
  GroupPack& gp = *pack_ctx(vctx);

  const size_t Nparams = gp.x.size();
  for (size_t D = 0; D < Nparams; D++) {
    gp.x[D] = x[D];
  }

  auto& pattern_here = innerv ? gp.pattern_grad_inner : gp.pattern_grad;
  auto& work_here = innerv ? gp.work_inner_grad : gp.work_grad;

  *out_f += gp.fun.Forward(0, gp.x)[0];
  gp.fun.sparse_jac_rev(
    gp.x,
    pattern_here,
    gp.unused_pattern,
    "cppad",
    work_here);

  const size_t NoutGrad = pattern_here.nnz();
  const auto& cols = pattern_here.col();
  const auto& vals = pattern_here.val();
  for (size_t D = 0; D < NoutGrad; D++) {
    out_grad[cols[D]] += vals[D];
  }

  return 0;
}

int eval_hess(
  void* vctx,
  const double* x,
  const bool* inner,
  double* out_f,
  double* out_grad,
  double* out_hes,
  int* map) {

  const bool innerv = *inner;
  GroupPack& gp = *pack_ctx(vctx);

  const size_t Nparams = gp.x.size();
  for (size_t D = 0; D < Nparams; D++) {
    gp.x[D] = x[D];
  }

  *out_f += gp.fun.Forward(0, gp.x)[0];

  auto& pattern_here_grad = innerv ? gp.pattern_grad_inner : gp.pattern_grad;
  auto& work_here_grad = innerv ? gp.work_inner_grad : gp.work_grad;

  auto& pattern_here_hes = innerv ? gp.pattern_hessian_inner : gp.pattern_hessian;
  auto& work_here_hes = innerv ? gp.work_inner_hess : gp.work_hess;

  gp.fun.sparse_jac_rev(
    gp.x,
    pattern_here_grad,
    gp.unused_pattern,
    "cppad",
    work_here_grad);

  gp.fun.sparse_hes(
    gp.x,
    gp.w,
    pattern_here_hes,
    gp.unused_pattern,
    "cppad.symmetric",
    work_here_hes);

  const size_t NoutGrad = pattern_here_grad.nnz();
  const auto& cols = pattern_here_grad.col();
  const auto& vals = pattern_here_grad.val();

  for (size_t D = 0; D < NoutGrad; D++) {
    out_grad[cols[D]] += vals[D];
  }

  const size_t n_hes = pattern_here_hes.nnz();
  const auto& vals_hes = pattern_here_hes.val();

  for (size_t Di = 0; Di < n_hes; Di++) {
    out_hes[map[Di]] = vals_hes[Di];
  }
  return 0;
}

int get_sparse_sizes(
  void* vctx,
  int* n_inner,
  int* n_outer,
  int* nnz_grad_inner,
  int* nnz_grad_outer,
  int* nnz_hes_inner,
  int* nnz_hes_outer) {

  GroupPack& gp = *pack_ctx(vctx);
  return get_sizes(
    gp, n_inner, n_outer, nnz_grad_inner, nnz_grad_outer, nnz_hes_inner, nnz_hes_outer);
}

int get_sparse_pattern(
  void* vctx,
  int* pattern_grad_inner,
  int* pattern_grad_outer,
  int* pattern_hes_inner_row,
  int* pattern_hes_inner_col,
  int* pattern_hes_outer_row,
  int* pattern_hes_outer_col) {

  GroupPack& gp = *pack_ctx(vctx);
  return get_pattern(
    gp,
    pattern_grad_inner,
    pattern_grad_outer,
    pattern_hes_inner_row,
    pattern_hes_inner_col,
    pattern_hes_outer_row,
    pattern_hes_outer_col);
}
