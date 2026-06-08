#ifndef ADLAPLACE_EVAL_TRACE_HINV_T_HPP
#define ADLAPLACE_EVAL_TRACE_HINV_T_HPP

#include <cstddef>

const char* trace_hinv_t_strerror(int rc);

int eval_trace_hinv_t(
  void* vctx,
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
  double* out_trace);

#endif
