#ifdef UNDEF

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/adpack.hpp"


inline CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
build_pattern_from_R(
  const size_t N,
  const Rcpp::IntegerVector& index,
  Rcpp::Nullable<Rcpp::IntegerVector> column = R_NilValue
){
#ifdef DEBUG
  // index: 0-based, strictly increasing, in-bounds
  if (index.size() > 0) {
    int prev = index[0];
    if (prev < 0 || (size_t)prev >= N)
      Rf_error("index out of bounds");

    for (int k = 1; k < index.size(); ++k) {
      int cur = index[k];
      if (cur <= prev)
        Rf_error("index must be strictly increasing (sorted, no duplicates)");
      if ((size_t)cur >= N)
        Rf_error("index out of bounds");
      prev = cur;
    }
  }
#endif

  const size_t K = (size_t) index.size();

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> pattern;

  if (column.isNull()) {
    // ---- Gradient pattern: 1 x N with nnz = K at (0, index[k])
    pattern.resize(1, N, K);
    for (size_t k = 0; k < K; ++k)
      pattern.set(k, 0, (size_t) index[(int)k]);
  } else {
    // ---- Hessian pattern: N x N with nnz = K at (index[k], column[k])
    Rcpp::IntegerVector col(column);

#ifdef DEBUG
    if (col.size() != index.size())
      Rf_error("column must have same length as index");
    for (int k = 0; k < col.size(); ++k) {
      int c = col[k];
      if (c < 0 || (size_t)c >= N)
        Rf_error("column out of bounds");
    }
#endif

    pattern.resize(N, N, K);
    for (size_t k = 0; k < K; ++k)
      pattern.set(k,
        (size_t) index[(int)k],
        (size_t) col[(int)k]
      );
  }

  return pattern;
}


inline Rcpp::List sparsityList(
        std::vector<PatternList> & pattern
        ) {

    const size_t Ngroups = pattern.size();
    Rcpp::List result(Ngroups);

    for(size_t D=0;D<Ngroups;D++) {

        result[D] = Rcpp::List::create(
            Rcpp::_["grad"] = Rcpp::wrap(pattern[D].grad.col()), 
            Rcpp::_["grad_inner"] = Rcpp::wrap(col_grad_inner.col()),

            Rcpp::_["row_hess"] = Rcpp::wrap(pattern[D].hessian_upper.row()),
            Rcpp::_["col_hess"] = Rcpp::wrap(pattern[D].hessian_upper.col()),

            Rcpp::_["row_hess_inner"] = Rcpp::wrap(pattern[D].hessian_inner_upper.row()),
            Rcpp::_["col_hess_inner"] = Rcpp::wrap(pattern[D].hessian_inner_upper.col())

            );
    }

    return(result);

}

#endif
