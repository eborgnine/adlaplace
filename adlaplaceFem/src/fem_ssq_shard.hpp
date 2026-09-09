#ifndef ADLAPLACEFEM_FEM_SSQ_SHARD_HPP
#define ADLAPLACEFEM_FEM_SSQ_SHARD_HPP

// ad_shard for random_fem_ssq_* that answers from closed-form derivatives
// instead of sweeping the AD tape. See fem_ssq_analytic.hpp for the formulas.
//
// EvalShard is final, so this derives from ad_shard directly and reuses the
// free helpers pack_sparsity_sizes / get_pattern for the pattern queries.
// Include only from the one translation unit that already included
// adlaplace/eval_impl.hpp.

#include "adlaplace/backend.hpp"

#include "fem_ssq_analytic.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

class FemSsqShard final : public ad_shard {
public:
  FemSsqShard(AdTape &&p, std::shared_ptr<const femssq::Data> data)
      : ad_shard(std::move(p), nullptr), data_(std::move(data)) {
    setup();
  }

  int f(const double *x, double *out_f) override {
    if (pack.x.size() == 0) {
      return 4;
    }
    gather(x);
    load_state(0);
    femssq::sym_quad_all(*data_, ws_.gamma, ws_.quad);
    *out_f += femssq::value_from_quad(*data_, ws_);
    return 0;
  }

  int f_grad(const double *x, bool inner, double *out_f,
             double *out_grad) override {
    if (pack.x.size() == 0) {
      return 4;
    }
    gather(x);
    load_state(1);
    femssq::sym_quad_all(*data_, ws_.gamma, ws_.quad);
    femssq::sym_matvec(*data_, ws_.Q_x, ws_.gamma, ws_.Qg);
    *out_f += femssq::value_from_quad(*data_, ws_);
    scatter_grad(inner, out_grad);
    return 0;
  }

  int f_grad_hess(const double *x, bool inner, double *out_f, double *out_grad,
                  double *out_hes, int *map) override {
    if (pack.x.size() == 0) {
      return 4;
    }
    const femssq::Data &d = *data_;
    gather(x);
    load_state(2);
    femssq::sym_quad_all(d, ws_.gamma, ws_.quad);
    femssq::sym_matvec(d, ws_.Q_x, ws_.gamma, ws_.Qg);
    for (std::size_t k = 0; k < femssq::kNumTheta; ++k) {
      femssq::sym_matvec(d, ws_.dQ_x[k], ws_.gamma, ws_.dQg[k]);
    }

    *out_f += femssq::value_from_quad(d, ws_);
    scatter_grad(inner, out_grad);

    const femssq::PatternRcv &pat =
        inner ? pack.pattern_hessian_inner : pack.pattern_hessian;
    const femssq::HesSlots &slots = inner ? hes_inner_ : hes_outer_;
    const std::size_t nnz = pat.nnz();
    for (std::size_t D = 0; D < nnz; ++D) {
      double v = 0.0;
      switch (slots.kind[D]) {
      case femssq::HesSlots::kGammaGamma:
        v = -ws_.Q_x[static_cast<std::size_t>(slots.a[D])];
        break;
      case femssq::HesSlots::kGammaTheta:
        v = -ws_.dQg[static_cast<std::size_t>(slots.b[D])]
                    [static_cast<std::size_t>(slots.a[D])];
        break;
      case femssq::HesSlots::kThetaTheta:
        v = femssq::hes_theta(d, ws_, static_cast<std::size_t>(slots.a[D]),
                              static_cast<std::size_t>(slots.b[D]));
        break;
      default:
        break;
      }
      // Assignment, not accumulation, matching EvalShard::f_grad_hess. Slots
      // with no analytic contribution must therefore be written as zeros.
      out_hes[map[D]] = v;
    }
    return 0;
  }

  int get_sizes(int *n_inner, int *n_outer, int *n_beta, int *n_theta,
                int *nnz_grad_inner, int *nnz_grad_outer, int *nnz_hes_inner,
                int *nnz_hes_outer) override {
    return pack_sparsity_sizes(pack, n_inner, n_outer, n_beta, n_theta,
                               nnz_grad_inner, nnz_grad_outer, nnz_hes_inner,
                               nnz_hes_outer);
  }

  int get_sparse_pattern(int *pattern_grad_inner, int *pattern_grad_outer,
                         int *pattern_hes_inner_row, int *pattern_hes_inner_col,
                         int *pattern_hes_outer_row,
                         int *pattern_hes_outer_col) override {
    return get_pattern(pack, pattern_grad_inner, pattern_grad_outer,
                       pattern_hes_inner_row, pattern_hes_inner_col,
                       pattern_hes_outer_row, pattern_hes_outer_col);
  }

  // No Forward/Reverse sweeps happen here, so there is no Taylor workspace
  // and no capacity_order to reserve. Do not drain pack.fun: a taped analytic
  // pack still has Domain() > 0, and capacity_order from the wrong DSO
  // trips CppAD thread_alloc.
  int assign_memory() override {
    return pack.x.size() == 0 ? 2 : 0;
  }

  void release_eval_buffers() override {
    pack.trace.direction.clear();
    pack.trace.direction_zeros.clear();
    pack.trace.wthree.clear();
  }

  int trace_hinv_t(const double *x, const int *LinvPt_p, const int *LinvPt_i,
                   const double *LinvPt_x, std::size_t LinvPt_ncol,
                   std::size_t LinvPt_p_len, std::size_t LinvPt_i_len,
                   std::size_t LinvPt_x_len, const int *LinvPtColumns_p,
                   const int *LinvPtColumns_i, std::size_t LinvPtColumns_p_len,
                   std::size_t LinvPtColumns_i_len,
                   double *out_trace) override {

    if (x == nullptr || out_trace == nullptr) {
      return 1;
    }
    if (LinvPt_p == nullptr || LinvPt_i == nullptr || LinvPt_x == nullptr) {
      return 1;
    }
    if (LinvPtColumns_p == nullptr || LinvPtColumns_i == nullptr) {
      return 1;
    }

    AdTape &gp = pack;
    const std::size_t ist = gp.shard_index;
    if (ist + 1 >= LinvPtColumns_p_len) {
      return 5;
    }
    if (LinvPt_p_len < LinvPt_ncol + 1) {
      return 6;
    }
    const std::size_t n_params = gp.x.size();
    if (n_params == 0) {
      return 2;
    }
    const bool compacted =
        gp.n_global > n_params && gp.tape_to_global.size() == n_params;

    const int col_start = LinvPtColumns_p[ist];
    const int col_end = LinvPtColumns_p[ist + 1];
    if (col_start < 0 || col_end < col_start) {
      return 7;
    }
    if (static_cast<std::size_t>(col_end) > LinvPtColumns_i_len) {
      return 8;
    }

    const femssq::Data &d = *data_;
    gather(x);
    load_state(1);

    double total[femssq::kNumTheta] = {0.0, 0.0};

    for (int dp = col_start; dp < col_end; ++dp) {
      const int dcol = LinvPtColumns_i[dp];
      if (dcol < 0 || static_cast<std::size_t>(dcol) >= LinvPt_ncol) {
        return 9;
      }
      if (static_cast<std::size_t>(dcol) + 1 >= LinvPt_p_len) {
        return 10;
      }
      const int entry_start = LinvPt_p[dcol];
      const int entry_end = LinvPt_p[dcol + 1];
      if (entry_start < 0 || entry_end < entry_start) {
        return 11;
      }
      if (static_cast<std::size_t>(entry_end) > LinvPt_i_len) {
        return 12;
      }
      if (static_cast<std::size_t>(entry_end) > LinvPt_x_len) {
        return 13;
      }

      // Stamp the support of this direction. Stale v_dense values are never
      // read because a row only counts when its stamp is current.
      ++ws_.stamp_tag;
      ws_.touched.clear();
      for (int k = entry_start; k < entry_end; ++k) {
        const int row = LinvPt_i[k];
        std::size_t tape_idx = 0;
        if (compacted) {
          if (row < 0 ||
              static_cast<std::size_t>(row) >= gp.gamma_row_to_tape.size()) {
            return 14;
          }
          tape_idx = gp.gamma_row_to_tape[static_cast<std::size_t>(row)];
          if (tape_idx >= n_params) {
            continue;
          }
        } else {
          if (row < 0 ||
              static_cast<std::size_t>(row) >= gp.pattern_grad_inner.nc()) {
            return 14;
          }
          tape_idx = gp.n_beta + static_cast<std::size_t>(row);
          if (tape_idx >= n_params) {
            return 15;
          }
        }
        const int local = t2g_[tape_idx];
        if (local < 0) {
          continue;
        }
        const std::size_t l = static_cast<std::size_t>(local);
        ws_.v_dense[l] = LinvPt_x[k];
        ws_.stamp[l] = ws_.stamp_tag;
        ws_.touched.push_back(local);
      }

      // v' (dQ/dtheta_t) v restricted to supp(v). Every upper-triangle pair
      // (r, c) with r <= c lives in CSC column c, so iterating the touched
      // columns and testing the row stamp visits each pair exactly once.
      for (int col_local : ws_.touched) {
        const std::size_t col = static_cast<std::size_t>(col_local);
        const double vc = ws_.v_dense[col];
        for (int pos = d.Q_p[col]; pos < d.Q_p[col + 1]; ++pos) {
          const std::size_t kk = static_cast<std::size_t>(pos);
          const std::size_t row = static_cast<std::size_t>(d.Q_i[kk]);
          if (ws_.stamp[row] != ws_.stamp_tag) {
            continue;
          }
          const double fac = d.w[kk] * vc * ws_.v_dense[row];
          for (std::size_t t = 0; t < femssq::kNumTheta; ++t) {
            total[t] += fac * ws_.dQ_x[t][kk];
          }
        }
      }
    }

    // f is quadratic in gamma, so the third derivative is zero on every gamma
    // coordinate; only the two thetas pick up a contribution.
    for (std::size_t t = 0; t < femssq::kNumTheta; ++t) {
      const std::size_t tape_idx = d.t_tape[t];
      const std::size_t g =
          compacted ? gp.tape_to_global[tape_idx] : tape_idx;
      out_trace[g] += -0.5 * total[t];
    }
    return 0;
  }

  ad_shard *clone() const override {
    return new FemSsqShard(clone_group_pack(pack), data_);
  }

private:
  std::shared_ptr<const femssq::Data> data_;
  std::vector<int> t2g_;
  femssq::GradSlots grad_outer_;
  femssq::GradSlots grad_inner_;
  femssq::HesSlots hes_outer_;
  femssq::HesSlots hes_inner_;
  femssq::Workspace ws_;

  void setup() {
    const femssq::Data &d = *data_;
    // Classify the patterns the pack already holds, so the nnz ordering is
    // correct by construction rather than by assuming how the analytic
    // pattern was enumerated.
    t2g_ = femssq::tape_to_gamma_row(d, pack.x.size());
    const auto qmap = femssq::build_q_index(d);
    grad_outer_ = femssq::build_grad_slots(d, pack.pattern_grad, t2g_);
    grad_inner_ = femssq::build_grad_slots(d, pack.pattern_grad_inner, t2g_);
    hes_outer_ = femssq::build_hes_slots(d, pack.pattern_hessian, t2g_, qmap);
    hes_inner_ =
        femssq::build_hes_slots(d, pack.pattern_hessian_inner, t2g_, qmap);
    ws_.init(d);
  }

  // Same gather as EvalShard: tape order, honouring compact_tape.
  void gather(const double *x) {
    AdTape &gp = pack;
    const std::size_t np = gp.x.size();
    if (gp.tape_to_global.size() == np) {
      for (std::size_t i = 0; i < np; ++i) {
        gp.x[i] = x[gp.tape_to_global[i]];
      }
    } else {
      for (std::size_t i = 0; i < np; ++i) {
        gp.x[i] = x[i];
      }
    }
  }

  void load_state(int order) {
    const femssq::Data &d = *data_;
    ws_.gamma.assign(d.n, 0.0);
    for (std::size_t r = 0; r < d.n; ++r) {
      ws_.gamma[r] = pack.x[d.gidx[r]];
    }
    double theta[femssq::kNumTheta];
    for (std::size_t k = 0; k < femssq::kNumTheta; ++k) {
      theta[k] = pack.x[d.t_tape[k]];
    }
    femssq::update_theta(d, ws_, theta, order);
  }

  void scatter_grad(bool inner, double *out_grad) {
    const femssq::Data &d = *data_;
    const femssq::PatternRcv &pat =
        inner ? pack.pattern_grad_inner : pack.pattern_grad;
    const femssq::GradSlots &slots = inner ? grad_inner_ : grad_outer_;
    const std::vector<std::size_t> &cols_global =
        inner ? pack.grad_inner_cols_global : pack.grad_cols_global;
    const std::size_t nnz = pat.nnz();
    const auto &cols = pat.col();
    const bool use_global = cols_global.size() == nnz;
    for (std::size_t D = 0; D < nnz; ++D) {
      const int target = slots.target[D];
      double v = 0.0;
      if (target >= 0) {
        v = -ws_.Qg[static_cast<std::size_t>(target)];
      } else if (target != femssq::GradSlots::kIgnore) {
        v = femssq::grad_theta(d, ws_,
                               static_cast<std::size_t>(-1 - target));
      }
      out_grad[use_global ? cols_global[D] : cols[D]] += v;
    }
  }
};

#endif
