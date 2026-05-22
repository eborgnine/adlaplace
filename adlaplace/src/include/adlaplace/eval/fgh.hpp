#ifndef ADLAPLACE_EVAL_FGH_HPP
#define ADLAPLACE_EVAL_FGH_HPP

#include "adlaplace/api/adpack_handle.h"

int eval_f(void* vctx, const double* x, double* out_f);

int eval_grad(void* vctx, const double* x, const bool* inner, double* out_f, double* out_grad);

int eval_hess(
  void* vctx,
  const double* x,
  const bool* inner,
  double* out_f,
  double* out_grad,
  double* out_hes,
  int* map);

int get_hessian(
  void* vctx,
  const bool* inner,
  const int** p,
  size_t* p_len,
  const int** i,
  size_t* i_len,
  int** map,
  size_t* map_len);

int get_sparse_sizes(
  void* vctx,
  int* n_inner,
  int* n_outer,
  int* nnz_grad_inner,
  int* nnz_grad_outer,
  int* nnz_hes_inner,
  int* nnz_hes_outer);

int get_sparse_pattern(
  void* vctx,
  int* pattern_grad_inner,
  int* pattern_grad_outer,
  int* pattern_hes_inner_row,
  int* pattern_hes_inner_col,
  int* pattern_hes_outer_row,
  int* pattern_hes_outer_col);

#endif
