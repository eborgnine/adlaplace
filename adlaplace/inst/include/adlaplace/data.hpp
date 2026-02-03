#ifndef DATACONFIG_HPP
#define DATACONFIG_HPP

#include<Rcpp.h>
#include<vector>
#include<cstddef>


// ----- safe getters -----
inline bool get_bool(const Rcpp::List& cfg, const char* key, bool def=false) {
  return cfg.containsElementNamed(key) ? Rcpp::as<bool>(cfg[key]) : def;
}
// helper to fetch int from R list with default
inline int get_int(const Rcpp::List& cfg, const char* key, int def = 1) {
  return cfg.containsElementNamed(key) ? Rcpp::as<int>(cfg[key]) : def;
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


// for vectors which might be integer or numeric

struct NumVecView {
  SEXP sexp = R_NilValue;
  bool is_int = false;

  NumVecView() {}

  explicit NumVecView(SEXP x) : sexp(x) {
    const int t = TYPEOF(x);
    if (t == INTSXP) {
      is_int = true;
    } else if (t == REALSXP) {
      is_int = false;
    } else {
      Rcpp::stop("Expected integer or numeric vector");
    }
  }

  inline size_t size() const {
    return static_cast<size_t>(XLENGTH(sexp));
  }

  inline double operator[](size_t i) const {
    return is_int
    ? static_cast<double>(INTEGER(sexp)[i])
    : REAL(sexp)[i];
  }
};


struct CscPattern {
  std::vector<std::size_t> i;   // row indices, length nnz
  std::vector<std::size_t> p;   // col pointers, length ncol+1
  std::size_t nrow = 0;
  std::size_t ncol = 0;

  CscPattern() = default;

  // Construct from Matrix::ngCMatrix / dgCMatrix / lgCMatrix (S4)
  explicit CscPattern(const Rcpp::S4& sm) {
    Rcpp::IntegerVector Ii = sm.slot("i");
    Rcpp::IntegerVector Pp = sm.slot("p");
    Rcpp::IntegerVector Dim = sm.slot("Dim");

    nrow = static_cast<std::size_t>(Dim[0]);
    ncol = static_cast<std::size_t>(Dim[1]);

    // copy into std::vector<size_t> (thread-safe thereafter)
    i.resize(static_cast<std::size_t>(Ii.size()));
    for (R_xlen_t k = 0; k < Ii.size(); ++k) {
      int v = Ii[k];
      i[static_cast<std::size_t>(k)] = static_cast<std::size_t>(v);
    }

    p.resize(static_cast<std::size_t>(Pp.size()));
    for (R_xlen_t k = 0; k < Pp.size(); ++k) {
      int v = Pp[k];
      p[static_cast<std::size_t>(k)] = static_cast<std::size_t>(v);
    }
  }
};

struct MatchGroup {
  CscPattern grad, grad_inner, hessian, hessian_inner;

  MatchGroup() = default;

  explicit MatchGroup(const Rcpp::List& obj) {

    grad        = CscPattern(Rcpp::S4(obj["grad"]));
    grad_inner  = CscPattern(Rcpp::S4(obj["grad_inner"]));
    hessian     = CscPattern(Rcpp::S4(obj["hessian"]));
    hessian_inner = CscPattern(Rcpp::S4(obj["hessian_inner"]));
  }
};

struct Data {
  DgCView A;
  DgCView X;

  NumVecView Qdiag, y;

  DgCView map;
  DgCView elgm_matrix;

  size_t Nmap, Nbeta, Ngamma, Ny;

  explicit Data(const Rcpp::List& data)
  : A(         DgCView(Rcpp::S4(data["ATp"])) )
  , X(         DgCView(Rcpp::S4(data["XTp"])) )
  , Qdiag(     data["Qdiag"] )
  , y(         data["y"] )
  , map(       DgCView(Rcpp::S4(data["map"])) )
  {

    if(data.containsElementNamed("elgm_matrix")) {
      elgm_matrix = DgCView(Rcpp::S4(data["elgm_matrix"]));
    } 
    Nmap = map.ncol();
    Nbeta   = static_cast<std::size_t>(X.nrow());
    Ngamma  = static_cast<std::size_t>(A.nrow());
    Ny    = static_cast<std::size_t>(y.size());   // == A.ncol()
    if(Ny != A.ncol() && A.nrow() != 0) {
      Rcpp::Rcout << "Ny " << Ny << " columns of A " << A.ncol() << "\n";
    }
    if(Ny != X.ncol() && X.nrow() != 0) {
      Rcpp::Rcout << "lengh y " << Ny << " columns of X " << X.ncol() << "\n";
    }
  }
};

// ----- config bundle -----
struct Config {
  bool verbose;
  bool transform_theta;
  int num_threads=1;

  Rcpp::NumericVector beta, gamma, theta;

  size_t beta_begin, beta_end, Nbeta, 
  gamma_begin, gamma_end, Ngamma, 
  theta_begin, theta_end, Ntheta, 
  Nparams, Ngroups;

  std::vector<size_t> Sgamma;
  std::vector<double> params;

  DgCView groups, hessian, hessian_inner;

  MatchGroup match;


  explicit Config(const Rcpp::List& cfg)
  : verbose(get_bool(cfg, "verbose", false)),
  transform_theta(get_bool(cfg, "transform_theta", false)),
  num_threads(get_int(cfg, "num_threads", 1)),
  beta(cfg["beta"]),
  gamma(cfg["gamma"]),
  theta(cfg["theta"])
  {

    beta_begin = 0;
    Nbeta = beta.size();
    beta_end = Nbeta;

    gamma_begin = beta_end;
    Ngamma = gamma.size();
    gamma_end = gamma_begin + Ngamma;

    theta_begin = gamma_end;
    Ntheta = theta.size();
    theta_end = theta_begin + Ntheta;

    Nparams = Nbeta + Ntheta + Ngamma;

    Sgamma = std::vector<size_t>(Ngamma);
    params = std::vector<double>(Nparams);

    for(size_t D=0;D<Nbeta;D++) {
      params[D] = beta[D];
    }
    for(size_t D=0,Dgamma=gamma_begin;D<Ngamma;D++,Dgamma++) {
      Sgamma[D] = Dgamma;
      params[Dgamma] = gamma[D];
    }
    for(size_t D=0,Dtheta = theta_begin;D<Ntheta;D++,Dtheta++) {
      params[Dtheta] = theta[D];
    } 

    Ngroups = 0;
    if(cfg.containsElementNamed("groups")) {
      groups = DgCView(Rcpp::S4(cfg["groups"]));
      Ngroups = groups.ncol();
    }

    if(cfg.containsElementNamed("sparsity")) {

      Rcpp::List sparsity = cfg["sparsity"];

      if(sparsity.containsElementNamed("hessian")) {
        hessian = DgCView(Rcpp::S4(sparsity["hessian"]));
      } 
      if(sparsity.containsElementNamed("hessian_inner")) {
        hessian_inner = DgCView(Rcpp::S4(sparsity["hessian_inner"]));
      }     

      if(sparsity.containsElementNamed("match")) {
        Rcpp::List matchList = sparsity["match"];
        match = MatchGroup(matchList);
      } 
    }
  }
};


#endif
