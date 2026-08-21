#ifndef ADMVN_PMVN_ATOMIC_HPP
#define ADMVN_PMVN_ATOMIC_HPP

#include "mvn_tape.hpp"
#include "mvn_sparsity.hpp"

#include <cppad/cppad.hpp>
#include <string>
#include <vector>

namespace admvn {
namespace detail {

// Defined in pmvn_atomic.cpp (one TU). Do not put thread_local in this
// header: Clang/MinGW emit duplicate TLS-init symbols per TU that includes
// an inline thread_local definition.
MvnTape* pmvn_atomic_tape_from_id(size_t call_id);
void set_pmvn_atomic_tape(size_t call_id, MvnTape* tape);
void set_pmvn_atomic_tapes(MvnTape* p1, MvnTape* p2);

class atomic_pmvn : public CppAD::atomic_four<double> {
public:
  explicit atomic_pmvn(const std::string& name)
    : CppAD::atomic_four<double>(name) {}

  bool for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y) override {
    (void)call_id;
    type_y.resize(1);
    type_y[0] = type_x[0];
    for (size_t i = 1; i < type_x.size(); ++i) {
      if (type_x[i] == CppAD::variable_enum) {
        type_y[0] = CppAD::variable_enum;
      }
    }
    return true;
  }

  bool rev_depend(
    size_t call_id,
    const CppAD::vector<bool>& ident_zero_x,
    CppAD::vector<bool>& depend_x,
    const CppAD::vector<bool>& depend_y) override {
    (void)call_id;
  (void)ident_zero_x;
    depend_x.resize(ident_zero_x.size());
    for (size_t i = 0; i < ident_zero_x.size(); ++i) {
      depend_x[i] = depend_y[0];
    }
    return true;
  }

  // Dense Jacobian sparsity: every selected domain input affects Phi.
  // Without this, for_jac_sparsity can omit parameters that only enter
  // through the orthant CDF (SUN L*, a, b, c), yielding exact zero grads.
  bool jac_sparsity(
    size_t call_id,
    bool dependency,
    const CppAD::vector<bool>& ident_zero_x,
    const CppAD::vector<bool>& select_x,
    const CppAD::vector<bool>& select_y,
    CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
    (void)call_id;
    (void)dependency;
    (void)ident_zero_x;
    const size_t n = select_x.size();
    const size_t m = select_y.size();
    size_t nnz = 0;
    if (m > 0 && select_y[0]) {
      for (size_t j = 0; j < n; ++j) {
        nnz += select_x[j];
      }
    }
    pattern_out.resize(m, n, nnz);
    size_t k = 0;
    if (m > 0 && select_y[0]) {
      for (size_t j = 0; j < n; ++j) {
        if (select_x[j]) {
          pattern_out.set(k++, 0, j);
        }
      }
    }
    return true;
  }

  bool hes_sparsity(
    size_t call_id,
    const CppAD::vector<bool>& ident_zero_x,
    const CppAD::vector<bool>& select_x,
    const CppAD::vector<bool>& select_y,
    CppAD::sparse_rc<CppAD::vector<size_t>>& pattern_out) override {
    (void)call_id;
    (void)ident_zero_x;
    const size_t n = select_x.size();
    size_t nnz = 0;
    if (select_y.size() > 0 && select_y[0]) {
      for (size_t i = 0; i < n; ++i) {
        if (!select_x[i]) {
          continue;
        }
        for (size_t j = i; j < n; ++j) {
          nnz += select_x[j];
        }
      }
    }
    pattern_out.resize(n, n, nnz);
    size_t k = 0;
    if (select_y.size() > 0 && select_y[0]) {
      for (size_t i = 0; i < n; ++i) {
        if (!select_x[i]) {
          continue;
        }
        for (size_t j = i; j < n; ++j) {
          if (select_x[j]) {
            pattern_out.set(k++, i, j);
          }
        }
      }
    }
    return true;
  }

  bool forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty) override {
    (void)select_y;
    if (order_low > 0 || order_up > 0) {
      return false;
    }
    ty.resize(1);
    MvnTape* pmvn_atomic_tape = pmvn_atomic_tape_from_id(call_id);
    if (pmvn_atomic_tape == nullptr) {
      return false;
    }
    const std::size_t n = pmvn_atomic_tape->n;
    const std::size_t n_in = 3 * n + vech_size(n);
    if (tx.size() < n_in) {
      return false;
    }

    std::vector<double> upper(n);
    std::vector<double> mean(n);
    GenzPack genz;
    genz.scale.resize(n);
    genz.ch.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
      upper[i] = tx[i];
      mean[i] = tx[n + i];
      genz.scale[i] = tx[2 * n + i];
    }
    const std::size_t off_ch = 3 * n;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        genz.ch[i][j] = tx[off_ch + vech_index(i, j)];
      }
    }

    double* no_err = nullptr;
    ty[0] = eval_mvn_value_double(*pmvn_atomic_tape, upper, mean, genz, no_err);
    return true;
  }

  bool reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py) override {
    (void)call_id;
    (void)select_x;
    (void)ty;
    if (order_up > 0) {
      return false;
    }
    MvnTape* pmvn_atomic_tape = pmvn_atomic_tape_from_id(call_id);
    if (pmvn_atomic_tape == nullptr || py.size() < 1) {
      return false;
    }
    const std::size_t n = pmvn_atomic_tape->n;
    const std::size_t n_in = 3 * n + vech_size(n);
    if (tx.size() < n_in) {
      return false;
    }
    px.resize(n_in);
    for (std::size_t i = 0; i < n_in; ++i) {
      px[i] = 0.0;
    }

    std::vector<double> upper(n);
    std::vector<double> mean(n);
    GenzPack genz;
    genz.scale.resize(n);
    genz.ch.assign(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
      upper[i] = tx[i];
      mean[i] = tx[n + i];
      genz.scale[i] = tx[2 * n + i];
    }
    const std::size_t off_ch = 3 * n;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        genz.ch[i][j] = tx[off_ch + vech_index(i, j)];
      }
    }

    std::vector<std::vector<double>> sigma(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        double s = 0.0;
        for (std::size_t k = 0; k <= std::min(i, j); ++k) {
          s += genz.ch[i][k] * genz.ch[j][k];
        }
        sigma[i][j] = genz.scale[i] * s * genz.scale[j];
      }
    }

    // Analytic domain gradient w.r.t. the (upper, mean, scale, L) in tx;
    // fall back to AD-through-QMC when analytic is unavailable (non-value_only).
    std::vector<double> grad_domain =
      eval_mvn_domain_grad_auto(*pmvn_atomic_tape, upper, mean, genz);
    if (grad_domain.empty()) {
      grad_domain = eval_mvn_domain_grad(*pmvn_atomic_tape, upper, mean, sigma);
    }
    if (grad_domain.size() != n_in) {
      return false;
    }

    for (std::size_t i = 0; i < n_in; ++i) {
      px[i] += py[0] * grad_domain[i];
    }
    return true;
  }
};

inline atomic_pmvn& pmvn_atomic_instance() {
  static atomic_pmvn op("admvn_pmvn");
  return op;
}

inline void init_pmvn_atomic() {
  (void)pmvn_atomic_instance();
}

inline CppAD::AD<double> pmvn_ad(
  size_t call_id,
  MvnTape& tape,
  const std::vector<CppAD::AD<double>>& upper,
  const std::vector<CppAD::AD<double>>& mean,
  const std::vector<CppAD::AD<double>>& scale,
  const std::vector<CppAD::AD<double>>& vech_ch) {

  const std::size_t n = tape.n;
  const std::size_t n_in = 3 * n + vech_size(n);
  CppAD::vector<CppAD::AD<double>> ax(n_in);
  for (std::size_t i = 0; i < n; ++i) {
    ax[i] = upper[i];
    ax[n + i] = mean[i];
    ax[2 * n + i] = scale[i];
  }
  const std::size_t off_ch = 3 * n;
  for (std::size_t i = 0; i < vech_ch.size(); ++i) {
    ax[off_ch + i] = vech_ch[i];
  }
  CppAD::vector<CppAD::AD<double>> ay(1);
  pmvn_atomic_instance()(call_id, ax, ay);
  return ay[0];
}

}  // namespace detail
}  // namespace admvn

#endif
