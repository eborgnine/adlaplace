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

// API-style third-derivative contraction used by older logLik derivative code.
// Returns 0 on success and nonzero on invalid inputs.
inline int eval_trace_hinv_t(
  void* vctx,
  const int* i,
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
  if (vctx == nullptr || i == nullptr || x == nullptr || out_trace == nullptr) return 1;
  if (LinvPt_p == nullptr || LinvPt_i == nullptr || LinvPt_x == nullptr) return 1;
  if (LinvPtColumns_p == nullptr || LinvPtColumns_i == nullptr) return 1;

  auto* ctx = static_cast<BackendContext*>(vctx);
  if (ctx->adFun == nullptr) return 2;
  if (*i < 0) return 3;

  const std::size_t ist = static_cast<std::size_t>(*i);
  if (ist >= ctx->adFun->size()) return 4;
  if (ist + 1 >= LinvPtColumns_p_len) return 5;
  if (LinvPt_p_len < LinvPt_ncol + 1) return 6;

  const CscView LinvPt{
    LinvPt_p,
    LinvPt_i,
    LinvPt_x,
    LinvPt_ncol,
    LinvPt_p_len,
    LinvPt_i_len,
    LinvPt_x_len
  };
  const CscView LinvPtColumns{
    LinvPtColumns_p,
    LinvPtColumns_i,
    nullptr,
    ctx->Ngroups,
    LinvPtColumns_p_len,
    LinvPtColumns_i_len,
    0
  };

  GroupPack& gp = (*(ctx->adFun))[ist];
  const std::size_t n_params = gp.x.size();
  std::copy_n(x, n_params, gp.x.begin());

  const int col_start = LinvPtColumns.p[ist];
  const int col_end = LinvPtColumns.p[ist + 1];
  if (col_start < 0 || col_end < col_start) return 7;
  if (static_cast<std::size_t>(col_end) > LinvPtColumns.i_len) return 8;

  for (int dp = col_start; dp < col_end; ++dp) {
    const int dcol = LinvPtColumns.i[dp];
    if (dcol < 0 || static_cast<std::size_t>(dcol) >= LinvPt.ncol) return 9;
    if (static_cast<std::size_t>(dcol) + 1 >= LinvPt.p_len) return 10;

    std::fill(gp.direction.begin(), gp.direction.end(), 0.0);

    const int entry_start = LinvPt.p[dcol];
    const int entry_end = LinvPt.p[dcol + 1];
    if (entry_start < 0 || entry_end < entry_start) return 11;
    if (static_cast<std::size_t>(entry_end) > LinvPt.i_len) return 12;
    if (static_cast<std::size_t>(entry_end) > LinvPt.x_len) return 13;

    for (int k = entry_start; k < entry_end; ++k) {
      const int row = LinvPt.i[k];
      if (row < 0 || static_cast<std::size_t>(row) >= ctx->Ngamma) return 14;

      const std::size_t idx = ctx->Nbeta + static_cast<std::size_t>(row);
      if (idx >= gp.direction.size()) return 15;
      gp.direction[idx] = LinvPt.x[k];
    }

    gp.fun.Forward(0, gp.x);
    gp.fun.Forward(1, gp.direction);
    gp.fun.Forward(2, gp.direction_zeros);

    const CppAD::vector<double> dw = gp.fun.Reverse(3, gp.wthree);
    if (dw.size() < 3 * n_params) return 16;

    for (std::size_t d = 0; d < n_params; ++d) {
      out_trace[d] += dw[3 * d];
    }
  }

  return 0;
}

#endif
