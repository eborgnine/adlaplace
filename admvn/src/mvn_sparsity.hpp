#ifndef ADMVN_MVN_SPARSITY_HPP
#define ADMVN_MVN_SPARSITY_HPP

#include <cppad/cppad.hpp>
#include <vector>

namespace admvn {
namespace detail {

static const char* kJacColor = "cppad";
static const char* kHessColor = "cppad.symmetric";

inline CppAD::vector<double> to_cppad_vector(const std::vector<double>& v) {
  CppAD::vector<double> out(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    out[i] = v[i];
  }
  return out;
}

inline void setup_mvn_sparsity(
  CppAD::ADFun<double>& fun,
  const CppAD::vector<double>& x,
  const std::vector<std::size_t>& inner_subset,
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>& pattern_grad,
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>& pattern_grad_inner,
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>& pattern_hessian,
  CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>& pattern_hessian_inner,
  CppAD::sparse_jac_work& work_grad,
  CppAD::sparse_jac_work& work_inner_grad,
  CppAD::sparse_hes_work& work_hess,
  CppAD::sparse_hes_work& work_inner_hess,
  CppAD::sparse_rc<CppAD::vector<size_t>>& unused_pattern,
  CppAD::vector<double>& w) {

  const size_t n_domain = x.size();
  w.resize(1);
  w[0] = 1.0;

  CppAD::sparse_rc<CppAD::vector<size_t>> grad;
  CppAD::sparse_rc<CppAD::vector<size_t>> grad_inner;
  CppAD::sparse_rc<CppAD::vector<size_t>> hessian;
  CppAD::sparse_rc<CppAD::vector<size_t>> hessian_upper;
  CppAD::sparse_rc<CppAD::vector<size_t>> hessian_inner;
  CppAD::sparse_rc<CppAD::vector<size_t>> hessian_inner_upper;

  CppAD::sparse_rc<CppAD::vector<size_t>> pattern_in;
  pattern_in.resize(n_domain, n_domain, n_domain);
  for (size_t j = 0; j < n_domain; ++j) {
    pattern_in.set(j, j, j);
  }

  CppAD::vectorBool select_domain(n_domain);
  CppAD::vectorBool select_range(1);
  for (size_t j = 0; j < n_domain; ++j) {
    select_domain[j] = true;
  }
  select_range[0] = true;

  const bool transpose = false;
  const bool dependency = false;
  const bool internal_bool = false;

  fun.for_jac_sparsity(pattern_in, transpose, dependency, internal_bool, grad);
  fun.for_hes_sparsity(select_domain, select_range, internal_bool, hessian);

  std::vector<unsigned char> insubset(n_domain, 0);
  for (size_t v : inner_subset) {
    if (v < n_domain) {
      insubset[v] = 1;
    }
  }

  {
    const auto& full_cols = grad.col();
    size_t k_inner = 0;
    for (size_t j : full_cols) {
      k_inner += insubset[j];
    }
    grad_inner.resize(1, n_domain, k_inner);
    size_t t = 0;
    for (size_t j : full_cols) {
      if (insubset[j]) {
        grad_inner.set(t++, 0, j);
      }
    }
  }

  {
    const auto& r = hessian.row();
    const auto& c = hessian.col();
    const size_t k = hessian.nnz();

    size_t k_upper = 0;
    size_t k_inner = 0;
    size_t k_inner_upper = 0;

    for (size_t idx = 0; idx < k; ++idx) {
      const size_t i = r[idx];
      const size_t j = c[idx];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];
      if (is_upper) {
        ++k_upper;
      }
      if (is_inner) {
        ++k_inner;
      }
      if (is_upper && is_inner) {
        ++k_inner_upper;
      }
    }

    hessian_upper.resize(n_domain, n_domain, k_upper);
    hessian_inner.resize(n_domain, n_domain, k_inner);
    hessian_inner_upper.resize(n_domain, n_domain, k_inner_upper);

    size_t tu = 0;
    size_t ti = 0;
    size_t tiu = 0;
    for (size_t idx = 0; idx < k; ++idx) {
      const size_t i = r[idx];
      const size_t j = c[idx];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];
      if (is_upper) {
        hessian_upper.set(tu++, i, j);
      }
      if (is_inner) {
        hessian_inner.set(ti++, i, j);
      }
      if (is_upper && is_inner) {
        hessian_inner_upper.set(tiu++, i, j);
      }
    }
  }

  fun.Forward(0, x);

  pattern_grad = CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>(grad);
  fun.sparse_jac_rev(x, pattern_grad, grad, kJacColor, work_grad);

  pattern_grad_inner =
    CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>(grad_inner);
  fun.sparse_jac_rev(x, pattern_grad_inner, grad_inner, kJacColor, work_inner_grad);

  pattern_hessian =
    CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>(hessian_upper);
  fun.sparse_hes(x, w, pattern_hessian, hessian, kHessColor, work_hess);

  pattern_hessian_inner =
    CppAD::sparse_rcv<CppAD::vector<size_t>, CppAD::vector<double>>(hessian_inner_upper);
  fun.sparse_hes(x, w, pattern_hessian_inner, hessian_inner, kHessColor, work_inner_hess);

  unused_pattern = hessian;
}

}  // namespace detail
}  // namespace admvn

#endif
