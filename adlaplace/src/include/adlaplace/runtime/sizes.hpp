#ifndef ADLAPLACE_SIZES_HPP
#define ADLAPLACE_SIZES_HPP

#include "adlaplace/api/backend.hpp"

int get_sizes(
  GroupPack& ad_pack,
  int* n_inner,
  int* n_outer,
  int* nnz_grad_inner,
  int* nnz_grad_outer,
  int* nnz_hes_inner,
  int* nnz_hes_outer);

int get_pattern(
  GroupPack& ad_pack,
  int* pattern_grad_inner,
  int* pattern_grad_outer,
  int* pattern_hes_inner_row,
  int* pattern_hes_inner_col,
  int* pattern_hes_outer_row,
  int* pattern_hes_outer_col);

#endif
