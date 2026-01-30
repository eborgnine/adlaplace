
#ifndef DGCVIEW_HPP
#define DGCVIEW_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

using PatternPair =
  std::pair<
    CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>,
    CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
  >;

using PatternPairs = std::vector<PatternPair>;

// ----- safe getters -----
inline bool get_bool(const Rcpp::List& cfg, const char* key, bool def=false) {
  return cfg.containsElementNamed(key) ? Rcpp::as<bool>(cfg[key]) : def;
}
// helper to fetch int from R list with default
inline int get_int(const Rcpp::List& cfg, const char* key, int def = 1) {
  return cfg.containsElementNamed(key) ? Rcpp::as<int>(cfg[key]) : def;
}

inline std::vector<double> get_numvec_copy(const Rcpp::List& cfg, const char* key) {
  if (!cfg.containsElementNamed(key)) return {};
  Rcpp::NumericVector v = cfg[key];
  return std::vector<double>(v.begin(), v.end());
}

// ----- a lightweight view of a dgCMatrix (no copies; just SEXP handles) -----
// Lightweight view for Matrix::*gCMatrix (column-compressed)
// Works for dgCMatrix (numeric x) and ngCMatrix (no x: pattern-only).
struct DgCView {
  Rcpp::IntegerVector i;    // row indices (0-based)
  Rcpp::IntegerVector p;    // column pointers (length ncol+1)
  Rcpp::NumericVector x;    // may be length 0 for ngCMatrix
  Rcpp::IntegerVector Dim;  // c(nrow, ncol)
  bool has_x;               // true if numeric/logical 'x' present

  DgCView()
    : i(), p(), x(), Dim(Rcpp::IntegerVector::create(0, 0)), has_x(false)
  {}

  explicit DgCView(const Rcpp::S4& obj)
  : i(obj.slot("i")),
  p(obj.slot("p")),
  Dim(obj.slot("Dim"))
  {
    // ngCMatrix has no 'x' slot; others do (dgCMatrix: numeric, lgCMatrix: logical)
    if (obj.inherits("ngCMatrix")) {
      x = Rcpp::NumericVector();   // empty
      has_x = false;
    } else {
      // coerce any existing x (numeric/logical) to NumericVector
      x = Rcpp::as<Rcpp::NumericVector>(obj.slot("x"));
      has_x = true;
    }
  }

  inline int nrow() const { return Dim[0]; }
  inline int ncol() const { return Dim[1]; }
  inline R_xlen_t nnz() const { return i.size(); } // structural nnz

  // Get numeric value for kth stored nonzero; returns 1 for pattern matrices.
  template <class T = double>
  inline T value(R_xlen_t k) const {
    return has_x ? static_cast<T>(x[k]) : static_cast<T>(1);
  }
};

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
        PatternPairs & pattern
        ) {

    const size_t Ngroups = pattern.size();
    Rcpp::List result(Ngroups);

    for(size_t D=0;D<Ngroups;D++) {

        const CppAD::vector<size_t>& col_grad = pattern[D].first.col();
        const CppAD::vector<size_t>& row_hess = pattern[D].second.row();
        const CppAD::vector<size_t>& col_hess = pattern[D].second.col();
        const size_t Kgrad = pattern[D].first.nnz();
        const size_t Khess = pattern[D].second.nnz();

        Rcpp::IntegerVector resultHereGrad(Kgrad);
        Rcpp::IntegerVector resultHereRow(Khess);
        Rcpp::IntegerVector resultHereCol(Khess);

        for (size_t k = 0; k < Kgrad; ++k) {
            resultHereGrad[k] = col_grad[k];
        }
        for (size_t k = 0; k < Khess; ++k) {
            resultHereRow[k] = row_hess[k];
            resultHereCol[k] = col_hess[k];
        }
        result[D] = Rcpp::List::create(
            Rcpp::_["grad"] = resultHereGrad, 
            Rcpp::_["row_hess"] = resultHereRow,
            Rcpp::_["col_hess"] = resultHereCol
            );
    }

    return(result);

}

#endif