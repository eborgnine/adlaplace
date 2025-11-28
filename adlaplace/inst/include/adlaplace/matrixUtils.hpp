#ifndef MATRIXUTILS_HPP
#define MATRIXUTILS_HPP

#include<Rcpp.h>
#include <RcppEigen.h>


// ----- safe getters -----
inline bool get_bool(const Rcpp::List& cfg, const char* key, bool def=false) {
  return cfg.containsElementNamed(key) ? Rcpp::as<bool>(cfg[key]) : def;
}
// helper to fetch int from R list with default
inline int get_int(const Rcpp::List& cfg, const char* key, int def = 1) {
  return cfg.containsElementNamed(key) ? Rcpp::as<int>(cfg[key]) : def;
}
inline double get_double(const Rcpp::List& cfg, const char* key, double def = 0.0) {
  return cfg.containsElementNamed(key) ? Rcpp::as<double>(cfg[key]) : def;
}

inline std::vector<double> get_numvec_copy(const Rcpp::List& cfg, const char* key) {
  if (!cfg.containsElementNamed(key)) return {};
  Rcpp::NumericVector v = cfg[key];
  return std::vector<double>(v.begin(), v.end());
}
inline std::vector<int> get_intvec_copy(const Rcpp::List& cfg, const char* key) {
  if (!cfg.containsElementNamed(key)) return {};
  Rcpp::IntegerVector v = cfg[key];
  return std::vector<int>(v.begin(), v.end());
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
// Lightweight view for Matrix::dsTMatrix (symmetric triplet), no uplo logic.
struct DgTView {
  Rcpp::IntegerVector i;    // row indices (0-based), as-is
  Rcpp::IntegerVector j;    // col indices (0-based), as-is
  Rcpp::NumericVector x;    // nonzeros, as-is
  Rcpp::IntegerVector Dim;  // c(nrow, ncol)

  explicit DgTView(const Rcpp::S4& obj)
  : i(obj.slot("i")),
  j(obj.slot("j")),
  x(obj.slot("x")),
  Dim(obj.slot("Dim"))
  {
    // (optional) assert type:
    // if (!obj.inherits("dsTMatrix")) Rcpp::stop("Expected dsTMatrix");
  }

  inline int nrow() const { return Dim[0]; }
  inline int ncol() const { return Dim[1]; }
  inline R_xlen_t nnz() const { return x.size(); }
};


inline Rcpp::IntegerVector compute_p_vector(
  const Rcpp::IntegerVector& j, 
  int ncol) {
    // Step 1: Count nonzeros per column
  std::vector<int> counts(ncol, 0);
  for (size_t k = 0; k < j.size(); ++k) {
        // j[k] is zero-based column index
    counts[j[k]]++;
  }

    // Step 2: Cumulative sum (CSC pointer)   
  Rcpp::IntegerVector p(ncol + 1);
  p[0] = 0;
  for (int c = 0; c < ncol; ++c) {
    p[c + 1] = p[c] + counts[c];
  }
//    return p;std::vector<int> p(ncol + 1, 0);
//    std::partial_sum(counts.begin(), counts.end(), p.begin() + 1);

  return p;
}



inline Rcpp::S4 make_TMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N)
{

  int Ni = N;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);

  Rcpp::LogicalVector toKeep = j >= 0;
  Rcpp::IntegerVector iR=i[toKeep], jR=j[toKeep];
  Rcpp::NumericVector xR=x[toKeep];


  Rcpp::S4 mat("dsTMatrix");
  mat.slot("i") = iR;
  mat.slot("j") = jR;
  mat.slot("x") = xR;
  mat.slot("Dim") = dims;
  mat.slot("uplo") =  Rcpp::wrap('L');
  return mat;



}

inline Rcpp::S4 make_CMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& p)
{

  int Ni = p.size()-1;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);

  Rcpp::S4 mat("dsCMatrix");
  mat.slot("i") = i;
  mat.slot("p") = p;
  mat.slot("x") = x;
  mat.slot("Dim") = dims;
  mat.slot("uplo") =  Rcpp::wrap('U');
  return mat;



}

inline Rcpp::S4 make_gCMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& p)
{

  int Ni = p.size()-1;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);

  Rcpp::S4 mat("dgCMatrix");
  mat.slot("i") = i;
  mat.slot("p") = p;
  mat.slot("x") = x;
  mat.slot("Dim") = dims;
  return mat;

}


inline Rcpp::RObject make_convert_gCmatrix(
  const std::vector<double>& x, 
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  const size_t N)
{
  const int nnz = i.size();

  std::vector<int> Iout;
  std::vector<int> Jout;
  std::vector<double> Xout;
  Iout.reserve(2 * nnz);
  Jout.reserve(2 * nnz);
  Xout.reserve(2 * nnz);
  for(int k=0;k<nnz;k++) {
    const int row = i[k];
    const int col = j[k];
    const double val = x[k]; 

    Iout.push_back(row);
    Jout.push_back(col);
    Xout.push_back(val);

    // add mirrored lower triangle if off-diagonal
    if (row != col) {
      Iout.push_back(col);
      Jout.push_back(row);
      Xout.push_back(val);
    }
  }
  const int nnzOut = Jout.size();
  std::vector<size_t> idx(nnzOut);
  std::iota(idx.begin(), idx.end(), 0);  

  std::stable_sort(idx.begin(), idx.end(),
    [&](size_t a, size_t b) { return Jout[a] < Jout[b]; });

  Rcpp::IntegerVector i1(nnzOut), j1(nnzOut);
  Rcpp::NumericVector x1(nnzOut);    
  for(int D=0;D<nnzOut;D++) {
    const int newD = idx[D];
    i1[D] = Iout[newD];
    j1[D] = Jout[newD];
    x1[D] = Xout[newD];
  }
  Rcpp::IntegerVector p1 = compute_p_vector(j1, N);
  Rcpp::S4 result = make_gCMatrix(x1, i1, p1);
  return(Rcpp::RObject(result));
}

// assumes you already have:
// inline Rcpp::S4 make_gCMatrix(
//   const Rcpp::NumericVector& x,
//   const Rcpp::IntegerVector& i,
//   const Rcpp::IntegerVector& p);

inline Rcpp::S4 eigen_to_dgC(const Eigen::SparseMatrix<double> &M) {
    using Eigen::Index;

    const Index nnz  = M.nonZeros();
    const Index ncol = M.cols();          // = outerSize() for ColMajor

    Rcpp::NumericVector x(nnz);
    Rcpp::IntegerVector i(nnz);
    Rcpp::IntegerVector p(ncol + 1);

    // copy from Eigen into R vectors
    std::copy(M.valuePtr(),
              M.valuePtr() + nnz,
              x.begin());

    std::copy(M.innerIndexPtr(),
              M.innerIndexPtr() + nnz,
              i.begin());

    std::copy(M.outerIndexPtr(),
              M.outerIndexPtr() + (ncol + 1),
              p.begin());

    // reuse your helper to build the S4 dgCMatrix
    return make_gCMatrix(x, i, p);
}

// Construct an N x N identity matrix as dgCMatrix layout.
inline DgCView identityMatrix(std::size_t N)
{
    const int n = static_cast<int>(N);

    const Rcpp::IntegerVector i = Rcpp::seq(0, n-1);       // row indices
    const Rcpp::IntegerVector p = Rcpp::seq(0, n);   // column pointers
    const Rcpp::NumericVector x = Rcpp::rep(1.0, n);  // all ones on the diagonal
    const Rcpp::IntegerVector Dim = Rcpp::rep(n, 2);

    Rcpp::S4 mat("dgCMatrix");
    mat.slot("i")   = i;
    mat.slot("p")   = p;
    mat.slot("x")   = x;
    mat.slot("Dim") = Dim;

    return DgCView(mat);
}

#endif

