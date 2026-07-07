#ifndef ADLAPLACE_EVAL_TRACE_HINV_T_EVAL_IMPL_HPP
#define ADLAPLACE_EVAL_TRACE_HINV_T_EVAL_IMPL_HPP

#include "adlaplace/eval/trace_hinv_t.hpp"

#include <algorithm>
#include <cppad/cppad.hpp>
#include <cstddef>
#include <vector>

#include "adlaplace/api/backend.hpp"

namespace {

struct CscView {
  const int* p;
  const int* i;
  const double* x;
  std::size_t ncol;
  std::size_t p_len;
  std::size_t i_len;
  std::size_t x_len;
};

}  // namespace

const char* trace_hinv_t_strerror(int rc) {
  switch (rc) {
    case 0: return "success";
    case 1: return "null input";
    case 2: return "empty parameter vector";
    case 5: return "LinvPtColumns p too short for shard_index";
    case 6: return "LinvPt p too short";
    case 7: return "invalid column range in LinvPtColumns";
    case 8: return "LinvPtColumns i index out of range";
    case 9: return "LinvPt column index out of range";
    case 10: return "LinvPt p too short for column";
    case 11: return "invalid LinvPt entry range";
    case 12: return "LinvPt i index out of range";
    case 13: return "LinvPt x index out of range";
    case 14: return "LinvPt row index out of range (gamma)";
    case 15: return "direction index out of range";
    case 16: return "Reverse(3) result too short";
    case 17: return "trace buffers not assigned";
    case 18: return "assign_memory API missing";
    default: return "unknown error";
  }
}

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
  double* out_trace) {

  if (vctx == nullptr || x == nullptr || out_trace == nullptr) return 1;
  if (LinvPt_p == nullptr || LinvPt_i == nullptr || LinvPt_x == nullptr) return 1;
  if (LinvPtColumns_p == nullptr || LinvPtColumns_i == nullptr) return 1;

  GroupPack& gp = *pack_ctx(vctx);
  const std::size_t ist = gp.shard_index;
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
    LinvPtColumns_p_len > 0 ? LinvPtColumns_p_len - 1 : 0,
    LinvPtColumns_p_len,
    LinvPtColumns_i_len,
    0
  };

  const std::size_t n_params = gp.fun.Domain();
  if (n_params == 0) return 2;

#ifdef DEBUG
  if (gp.direction.size() != n_params ||
      gp.direction_zeros.size() != n_params ||
      gp.wthree.size() != 3) {
    return 17;
  }
#endif

  int rc = 0;
  if (gp.x.size() != n_params) {
    rc = 2;
  } else {
    std::copy_n(x, n_params, gp.x.begin());

    const int col_start = LinvPtColumns.p[ist];
    const int col_end = LinvPtColumns.p[ist + 1];
    if (col_start < 0 || col_end < col_start) {
      rc = 7;
    } else if (static_cast<std::size_t>(col_end) > LinvPtColumns.i_len) {
      rc = 8;
    } else {
      for (int dp = col_start; dp < col_end && rc == 0; ++dp) {
        const int dcol = LinvPtColumns.i[dp];
        if (dcol < 0 || static_cast<std::size_t>(dcol) >= LinvPt.ncol) {
          rc = 9;
          break;
        }
        if (static_cast<std::size_t>(dcol) + 1 >= LinvPt.p_len) {
          rc = 10;
          break;
        }

        std::fill(gp.direction.begin(), gp.direction.end(), 0.0);

        const int entry_start = LinvPt.p[dcol];
        const int entry_end = LinvPt.p[dcol + 1];
        if (entry_start < 0 || entry_end < entry_start) {
          rc = 11;
          break;
        }
        if (static_cast<std::size_t>(entry_end) > LinvPt.i_len) {
          rc = 12;
          break;
        }
        if (static_cast<std::size_t>(entry_end) > LinvPt.x_len) {
          rc = 13;
          break;
        }

        for (int k = entry_start; k < entry_end; ++k) {
          const int row = LinvPt.i[k];
          if (row < 0 || static_cast<std::size_t>(row) >= gp.pattern_grad_inner.nc()) {
            rc = 14;
            break;
          }

          const std::size_t idx = gp.n_beta + static_cast<std::size_t>(row);
          if (idx >= gp.direction.size()) {
            rc = 15;
            break;
          }
          gp.direction[idx] = LinvPt.x[k];
        }
        if (rc != 0) {
          break;
        }

        gp.fun.Forward(0, gp.x);
        gp.fun.Forward(1, gp.direction);
        gp.fun.Forward(2, gp.direction_zeros);
        const CppAD::vector<double> dw = gp.fun.Reverse(3, gp.wthree);
        if (dw.size() < 3 * n_params) {
          rc = 16;
          break;
        }

        for (std::size_t d = 0; d < n_params; ++d) {
          out_trace[d] += dw[3 * d];
        }
      }
    }
  }

  gp.direction.clear();
  gp.direction_zeros.clear();
  gp.wthree.clear();
  return rc;
}

#endif
