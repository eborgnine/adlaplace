#pragma once

#ifndef LOGSPACE_HPOL_HPP
#define LOGSPACE_HPOL_HPP


#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <type_traits>

#include <utility>      // for std::pair


#include"lgamma.hpp"
#include"logspaceadd.hpp"
#include"matrixUtils.hpp"


// Constants
inline constexpr double LOGTWOPI      = 1.8378770664093454836;
inline constexpr double HALFLOGTWOPI  = 0.9189385332046727;


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


struct DataR {
  // existing
  DgTView QsansDiag;   // if this is actually dsTMatrix, use DgTView instead
  DgCView A;
  DgCView X;
  DgCView CC;

  Rcpp::NumericVector Qdiag, y;
  Rcpp::IntegerVector map;
  size_t Nq, Nbeta, Ngamma, Neta, Nstrata;

  // new (for dynamic payload planning)
  bool   has_offset = false;
  size_t Kmax_strata = 0;    // max NinStrata
  size_t Kx_col_max  = 0;    // max nnz in any X column
  size_t Ka_col_max  = 0;    // max nnz in any A column

  // convenience: dynamic payload sizes
  size_t Ndyn_dense_per_obs  = 0;
  size_t Ndyn_dense_total    = 0;
  size_t Ndyn_sparse_per_obs = 0;  // only valid if you bucket by fixed pattern
  size_t Ndyn_sparse_total   = 0;

  // (optionally) store ncol for X/A if your S4 wrapper supports it
  size_t X_ncol = 0;
  size_t A_ncol = 0;

  explicit DataR(const Rcpp::List& data)
  : QsansDiag( DgTView(Rcpp::S4(data["QsansDiag"])) )  
    , A(         DgCView(Rcpp::S4(data["ATp"])) )
    , X(         DgCView(Rcpp::S4(data["XTp"])) )
    , CC(        DgCView(Rcpp::S4(data["cc_matrixTp"])) )
    , Qdiag(     data["Qdiag"] )
    , y(         data["y"] )
    , map(       data["map"] )
  {
  // Basic sizes
    Nq      = static_cast<std::size_t>(QsansDiag.nnz());
    Nbeta   = static_cast<std::size_t>(X.nrow());
    Ngamma  = static_cast<std::size_t>(A.nrow());
    Neta    = static_cast<std::size_t>(X.ncol());   // == A.ncol()
    Nstrata = static_cast<std::size_t>(CC.ncol());

    has_offset = data.containsElementNamed("offset");
    X_ncol = static_cast<std::size_t>(X.ncol());
    A_ncol = static_cast<std::size_t>(A.ncol());

    // ---- Kmax over strata (max NinStrata) ----
    {
      for (size_t s = 0; s < Nstrata; ++s) {
        const size_t nin = static_cast<size_t>(CC.p[s+1] - CC.p[s]);
        if (nin > Kmax_strata) Kmax_strata = nin;
      }
    }

    // Max obs per stratum (column nnz of CC)
    for (std::size_t s = 0; s < Nstrata; ++s) {
      const std::size_t nin = static_cast<std::size_t>(CC.p[s+1] - CC.p[s]);
      if (nin > Kmax_strata) Kmax_strata = nin;
    }

    // Max nnz per column in X and A
    for (std::size_t c = 0; c < Neta; ++c) {
      const std::size_t nnzX = static_cast<std::size_t>(X.p[c+1] - X.p[c]);
      const std::size_t nnzA = static_cast<std::size_t>(A.p[c+1] - A.p[c]);
      if (nnzX > Kx_col_max) Kx_col_max = nnzX;
      if (nnzA > Ka_col_max) Ka_col_max = nnzA;
    }
   // Dynamic payload sizes
    Ndyn_dense_per_obs = 1 /*mask*/ + 1 /*y*/ + (has_offset ? 1 : 0) + Nbeta + Ngamma;
    Ndyn_dense_total   = Kmax_strata * Ndyn_dense_per_obs;

  }

};


// Thread-safe, R-free copies of your inputs
struct CscMatrix {
  // column-compressed sparse matrix (like Matrix::dgCMatrix / ngCMatrix)
  std::vector<int>    i;     // row indices (length = nnz)
  std::vector<int>    p;     // column pointers (length = ncol+1)
  std::vector<double> x;     // nonzero values (length = nnz); may be all 1s for pattern
  int nrow = 0, ncol = 0;

  inline size_t nnz() const { return i.size(); }
};

struct TripletMatrix {
  // (symmetric) triplet form like Matrix::dsTMatrix (we store as-is)
  std::vector<int>    i;     // row indices (0-based)
  std::vector<int>    j;     // col indices (0-based)
  std::vector<double> x;     // values
  int nrow = 0, ncol = 0;

  inline size_t nnz() const { return x.size(); }
};

struct Data {
  // Sparse matrices (owned, plain vectors)
  TripletMatrix QsansDiag;  // (off-diagonal part; you treat it like dsT)
  CscMatrix     A;          // rows = Ngamma, cols = Neta
  CscMatrix     X;          // rows = Nbeta,  cols = Neta
  CscMatrix     CC;         // selector (obs-by-strata)

  // Dense vectors (owned)
  std::vector<double> Qdiag;   // length Ngamma
  std::vector<double> y;       // length Neta (or #obs; per your code, it's indexed by CC.i)
  std::vector<int>    map;     // length Ngamma

  // Sizes (cached)
  size_t Nq       = 0;
  size_t Nbeta    = 0;
  size_t Ngamma   = 0;
  size_t Neta     = 0;
  size_t Nstrata  = 0;

  // Optional / derived
  bool   has_offset   = false;
  size_t Kmax_strata  = 0;   // max number of obs in any stratum (column nnz of CC)
  size_t Kx_col_max   = 0;   // max nnz in any X column
  size_t Ka_col_max   = 0;   // max nnz in any A column

  // (Optional) dynamic payload sizing you computed before
  size_t Ndyn_dense_per_obs  = 0;
  size_t Ndyn_dense_total    = 0;
  size_t Ndyn_sparse_per_obs = 0;
  size_t Ndyn_sparse_total   = 0;

  // ---------- helpers to copy from Rcpp S4 into plain structs ----------
  static CscMatrix copy_dgc(const Rcpp::S4& mat) {
    CscMatrix out;
    Rcpp::IntegerVector Iv = mat.slot("i");
    Rcpp::IntegerVector Pv = mat.slot("p");
    Rcpp::IntegerVector Dv = mat.slot("Dim");

    out.i.assign(Iv.begin(), Iv.end());
    out.p.assign(Pv.begin(), Pv.end());
    out.nrow = Dv[0];
    out.ncol = Dv[1];

    // Some Matrix classes (ngCMatrix) have no 'x'. Handle generically:
    if (mat.hasSlot("x")) {
      Rcpp::NumericVector Xv = mat.slot("x");
      out.x.assign(Xv.begin(), Xv.end());
    } else {
      // pattern-only matrix: treat all structural nonzeros as 1.0
      out.x.assign(out.i.size(), 1.0);
    }
    return out;
  }

  static TripletMatrix copy_dst(const Rcpp::S4& mat) {
    TripletMatrix out;
    Rcpp::IntegerVector Iv = mat.slot("i");
    Rcpp::IntegerVector Jv = mat.slot("j");
    Rcpp::NumericVector Xv = mat.slot("x");
    Rcpp::IntegerVector Dv = mat.slot("Dim");

    out.i.assign(Iv.begin(), Iv.end());
    out.j.assign(Jv.begin(), Jv.end());
    out.x.assign(Xv.begin(), Xv.end());
    out.nrow = Dv[0];
    out.ncol = Dv[1];
    return out;
  }

  // ---------- main constructor from the original R list ----------
  explicit Data(const Rcpp::List& data)
  : QsansDiag(copy_dst( Rcpp::S4(data["QsansDiag"]) ))
  , A(        copy_dgc( Rcpp::S4(data["ATp"])      ))
  , X(        copy_dgc( Rcpp::S4(data["XTp"])      ))
  , CC(       copy_dgc( Rcpp::S4(data["cc_matrixTp"]) ))
  {
    // Dense vectors
    {
      Rcpp::NumericVector rQdiag = data["Qdiag"];
      Qdiag.assign(rQdiag.begin(), rQdiag.end());

      Rcpp::NumericVector ry = data["y"];
      y.assign(ry.begin(), ry.end());

      Rcpp::IntegerVector rmap = data["map"];
      map.assign(rmap.begin(), rmap.end());
    }

    // Sizes (mirror your earlier logic)
    Nq       = QsansDiag.nnz();
    Nbeta    = static_cast<size_t>(X.nrow);
    Ngamma   = static_cast<size_t>(A.nrow);
    Neta     = static_cast<size_t>(X.ncol);   // == A.ncol (by construction)
    Nstrata  = static_cast<size_t>(CC.ncol);

    has_offset = data.containsElementNamed("offset");

    // ---- Kmax over strata (max nin per stratum = max column nnz of CC) ----
    Kmax_strata = 0;
    for (size_t s = 0; s < Nstrata; ++s) {
      const size_t nin = static_cast<size_t>(CC.p[s+1] - CC.p[s]);
      if (nin > Kmax_strata) Kmax_strata = nin;
    }

    // ---- Max nnz per column in X and A ----
    Kx_col_max = 0;
    Ka_col_max = 0;
    for (size_t c = 0; c < static_cast<size_t>(Neta); ++c) {
      const size_t nnzX = static_cast<size_t>(X.p[c+1] - X.p[c]);
      const size_t nnzA = static_cast<size_t>(A.p[c+1] - A.p[c]);
      if (nnzX > Kx_col_max) Kx_col_max = nnzX;
      if (nnzA > Ka_col_max) Ka_col_max = nnzA;
    }

    // ---- Optional: your payload sizing (unchanged) ----
    Ndyn_dense_per_obs = 1 /*mask*/ + 1 /*y*/ + (has_offset ? 1 : 0)
                       + Nbeta + Ngamma;
    Ndyn_dense_total   = Kmax_strata * Ndyn_dense_per_obs;

    Ndyn_sparse_per_obs = 1 /*mask*/ + 1 /*y*/ + (has_offset ? 1 : 0)
                        + (Kx_col_max + Ka_col_max);
    Ndyn_sparse_total   = Kmax_strata * Ndyn_sparse_per_obs;
  }

  // ---------- tiny convenience accessors (optional) ----------
  inline int X_col_nnz(int col) const { return pdiff(X.p, col); }
  inline int A_col_nnz(int col) const { return pdiff(A.p, col); }
  inline int CC_col_nnz(int col) const { return pdiff(CC.p, col); }

private:
  static inline int pdiff(const std::vector<int>& p, int c) {
    return p[c+1] - p[c];
  }
};

// ----- config bundle -----
struct Config {
  Rcpp::List list;             
  bool verbose;
  bool dirichlet;
  bool transform_theta;

  // optional numeric vectors (may be length 0)
  std::vector<double>beta;
  std::vector<double> theta;
  int num_threads=1;
  int maxDeriv=0;
  bool dense=0;
  Rcpp::List sparsity;
  double halfLogDetQ;

  explicit Config(const Rcpp::List& cfg)
  : list(cfg),
  verbose(get_bool(cfg, "verbose", false)),
  dirichlet(get_bool(cfg, "dirichlet", false)),
  transform_theta(get_bool(cfg, "transform_theta", false)),
  beta(get_numvec_copy(cfg, "beta")),
  theta(get_numvec_copy(cfg, "theta")),
  num_threads(get_int(cfg, "num_threads", 1)),          // <-- int, default 1
  maxDeriv(get_int(cfg, "maxDeriv", 0)),
  dense(get_bool(cfg, "dense", 0)),
  sparsity(cfg.containsElementNamed("sparsity") ? cfg["sparsity"] : Rcpp::List()),
  halfLogDetQ(get_double(cfg, "halfLogDetQ", 0))
  {}
};

template<class Type>
struct PackedParams {
  size_t Nbeta, Ngamma, Ntheta;
  bool dirichlet, transform_theta;
  CppAD::vector<Type> beta;     // size Nbeta
  CppAD::vector<Type> gamma;    // size Ngamma
  CppAD::vector<Type> theta;    // size Ntheta (on natural scale)
  CppAD::vector<Type> logTheta; // size Ntheta (log scale)
  Type logSqrtNu;
  Type oneOverSqrtNu;
  Type lgammaOneOverSqrtNu;
  size_t startGamma = 0;        // index into ad_params where gamma starts (if sourced there)
};


template<class Type>
inline PackedParams<Type>
unpack_params(const CppAD::vector<Type>& params,
              const Data& data,
              const Config& cfg)
{
  using CppAD::exp; using CppAD::log;

std::size_t Ntheta_base = 0u;

if (!data.map.empty()) {
    int max_val = data.map[0];
    for (std::size_t i = 1; i < data.map.size(); ++i) {
        if (data.map[i] > max_val) {
            max_val = data.map[i];
        }
    }
    Ntheta_base = static_cast<std::size_t>(max_val) + 1;
}


  PackedParams<Type> out;
  out.Nbeta = data.Nbeta;
  out.Ngamma = data.Ngamma;
  out.Ntheta = Ntheta_base + (cfg.dirichlet ? 1u : 0u);
  out.dirichlet = cfg.dirichlet;
  out.transform_theta = cfg.transform_theta;

  out.beta     = CppAD::vector<Type>(out.Nbeta);
  out.gamma    = CppAD::vector<Type>(out.Ngamma);
  out.theta    = CppAD::vector<Type>(out.Ntheta);
  out.logTheta = CppAD::vector<Type>(out.Ntheta);

  const std::size_t Nparams = params.size();

  if (Nparams == out.Ngamma) { // parameters are gamma only
    for (std::size_t Dparam=0; Dparam<out.Ngamma; ++Dparam) out.gamma[Dparam] = params[Dparam];

    if (static_cast<std::size_t>(cfg.beta.size()) != out.Nbeta)
      Rcpp::stop("beta has wrong length: expected %zu, got %d", out.Nbeta, cfg.beta.size());
    for (std::size_t Dparam=0; Dparam<out.Nbeta; ++Dparam) out.beta[Dparam] = Type(cfg.beta[Dparam]);

    if (static_cast<std::size_t>(cfg.theta.size()) != out.Ntheta)
      Rcpp::stop("theta has wrong length: expected %zu, got %d", out.Ntheta, cfg.theta.size());

    if (cfg.transform_theta) {
      for (std::size_t i=0; i<out.Ntheta; ++i) {
        out.logTheta[i] = Type(cfg.theta[i]);
        out.theta[i]    = exp(out.logTheta[i]);
      }
    } else {
      for (std::size_t i=0; i<out.Ntheta; ++i) {
        out.theta[i]    = Type(cfg.theta[i]);
        out.logTheta[i] = log(out.theta[i]);
      }
    }
  } else { // parameters is beta, gamma theta
    const std::size_t expected = out.Nbeta + out.Ngamma + out.Ntheta;
    if (Nparams != expected)
      Rcpp::stop("parameters has wrong size: expected %zu, got %zu", expected, Nparams);

    for (std::size_t i=0; i<out.Nbeta;  ++i) out.beta[i]  = params[i];
    out.startGamma = out.Nbeta;
    for (std::size_t i=0; i<out.Ngamma; ++i) out.gamma[i] = params[out.startGamma + i];

    const std::size_t startTheta = out.Nbeta + out.Ngamma;
    if (cfg.transform_theta) {
      for (std::size_t i=0; i<out.Ntheta; ++i) {
        out.logTheta[i] = params[startTheta + i];
        out.theta[i]    = exp(out.logTheta[i]);
      }
    } else {
      for (std::size_t i=0; i<out.Ntheta; ++i) {
        out.theta[i]    = params[startTheta + i];
        out.logTheta[i] = log(out.theta[i]);
      }
    }
  }

  if (cfg.dirichlet) {
    out.logSqrtNu           = out.logTheta[out.logTheta.size()-1] / Type(2);
    out.oneOverSqrtNu       = CppAD::exp(-out.logSqrtNu);
    out.lgammaOneOverSqrtNu = lgamma_ad(out.oneOverSqrtNu);
  } else {
    out.logSqrtNu = Type(0);
    out.oneOverSqrtNu = Type(0);
    out.lgammaOneOverSqrtNu = Type(0);
  }
  return out;
}

inline PackedParams<double>
unpack_params(const Rcpp::NumericVector& params,
              const Data& data,
              const Config& cfg)
{
  size_t Nparams = params.size();
  CppAD::vector<double> params2(Nparams);
  for(size_t D=0;D<Nparams;++D) {
    params2[D] = params[D];
  }

  PackedParams<double> result=unpack_params<double>(params2, data, cfg);
  return(result);
}


template<class Type>
CppAD::vector<Type>  objectiveFunctionInternal(
 const CppAD::vector<Type>& ad_params,  
 const Data& data,
 const Config& config
 );


#endif // LOGSPACE_HPOL_HPP