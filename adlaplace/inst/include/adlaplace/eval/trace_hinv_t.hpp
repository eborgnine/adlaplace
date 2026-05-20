#ifndef ADLAPLACE_CREATORS_THIRD_HPP
#define ADLAPLACE_CREATORS_THIRD_HPP

#include <algorithm>
#include <cppad/cppad.hpp>
#include <cstddef>
#include <vector>

#include "adlaplace/runtime/backend.hpp"

struct CscView {
  const int* p;
  const int* i;
  const double* x;
  std::size_t ncol;
  std::size_t p_len;
  std::size_t i_len;
  std::size_t x_len;
};

// Disabled until trace path is redesigned without BackendContext metadata.
inline int eval_trace_hinv_t(
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
  double* out_trace
) {
  (void)vctx;
  (void)x;
  (void)LinvPt_p;
  (void)LinvPt_i;
  (void)LinvPt_x;
  (void)LinvPt_ncol;
  (void)LinvPt_p_len;
  (void)LinvPt_i_len;
  (void)LinvPt_x_len;
  (void)LinvPtColumns_p;
  (void)LinvPtColumns_i;
  (void)LinvPtColumns_p_len;
  (void)LinvPtColumns_i_len;
  (void)out_trace;
  return 1;
}

#endif
