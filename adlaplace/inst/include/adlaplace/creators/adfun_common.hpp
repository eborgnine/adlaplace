#ifndef ADFUN_COMMON_HPP
#define ADFUN_COMMON_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/api/backend.hpp"

static const std::string JAC_COLOR  = "cppad";
static const std::string HESS_COLOR = "cppad.symmetric";

template <class SizeVec, class ValVec>
inline ValVec& rcv_val(CppAD::sparse_rcv<SizeVec, ValVec>& rcv) {
  return const_cast<ValVec&>(rcv.val());
}

inline void adpack_sparsity(
  const CPPAD_TESTVECTOR(double)& x,
  const std::vector<int>& subset,
  GroupPack& gp,
  const bool verbose = false) {
  const size_t Nparams = x.size();

  gp.w.resize(1);
  gp.w[0] = 1.0;

  gp.wthree.resize(3);
  gp.wthree[0] = 0.0;
  gp.wthree[1] = 0.0;
  gp.wthree[2] = 1.0;

  gp.x.resize(Nparams);
  gp.direction.resize(Nparams);
  gp.direction_zeros.resize(Nparams);

  std::fill(gp.direction.begin(), gp.direction.end(), 0.0);
  std::fill(gp.direction_zeros.begin(), gp.direction_zeros.end(), 0.0);

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad_inner;

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_upper;

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner_upper;

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> pattern_in;
  pattern_in.resize(Nparams, Nparams, Nparams);
  for (size_t j = 0; j < Nparams; ++j) {
    pattern_in.set(j, j, j);
  }

  CppAD::vectorBool select_domain(Nparams), select_range(1);
  for (size_t j = 0; j < Nparams; ++j) {
    select_domain[j] = true;
  }
  select_range[0] = true;

  const bool transpose     = false;
  const bool dependency    = false;
  const bool internal_bool = false;

  gp.fun.for_jac_sparsity(pattern_in, transpose, dependency, internal_bool, grad);
  gp.fun.for_hes_sparsity(select_domain, select_range, internal_bool, hessian);

  const auto& full_cols = grad.col();

  std::vector<unsigned char> insubset(Nparams, 0);
  for (size_t v : subset) {
    insubset[v] = 1;
  }

  {
    size_t K = 0;
    for (size_t j : full_cols) {
      K += insubset[j];
    }
    grad_inner.resize(1, Nparams, K);

    size_t t = 0;
    for (size_t j : full_cols) {
      if (insubset[j]) grad_inner.set(t++, 0, j);
    }
  }

  {
    const auto& r = hessian.row();
    const auto& c = hessian.col();
    const size_t K = hessian.nnz();

    size_t K_upper = 0, K_inner = 0, K_inner_upper = 0;

    for (size_t k = 0; k < K; ++k) {
      const size_t i = r[k], j = c[k];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper) ++K_upper;
      if (is_inner) ++K_inner;
      if (is_upper && is_inner) ++K_inner_upper;
    }

    hessian_upper.resize(Nparams, Nparams, K_upper);
    hessian_inner.resize(Nparams, Nparams, K_inner);
    hessian_inner_upper.resize(Nparams, Nparams, K_inner_upper);

    size_t tu = 0, ti = 0, tiu = 0;
    for (size_t k = 0; k < K; ++k) {
      const size_t i = r[k], j = c[k];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper) hessian_upper.set(tu++, i, j);
      if (is_inner) hessian_inner.set(ti++, i, j);
      if (is_upper && is_inner) hessian_inner_upper.set(tiu++, i, j);
    }
  }

  gp.fun.Forward(0, x);

  gp.pattern_grad = CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(grad);
  gp.fun.sparse_jac_rev(x, gp.pattern_grad, grad, JAC_COLOR, gp.work_grad);

  gp.pattern_grad_inner = CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(grad_inner);
  gp.fun.sparse_jac_rev(x, gp.pattern_grad_inner, grad_inner, JAC_COLOR, gp.work_inner_grad);

  gp.pattern_hessian = CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(hessian_upper);
  gp.fun.sparse_hes(x, gp.w, gp.pattern_hessian, hessian, HESS_COLOR, gp.work_hess);

  gp.pattern_hessian_inner = CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(hessian_inner_upper);
  gp.fun.sparse_hes(x, gp.w, gp.pattern_hessian_inner, hessian_inner, HESS_COLOR, gp.work_inner_hess);
}

#endif
