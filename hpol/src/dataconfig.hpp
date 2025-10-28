#ifndef DATA_HPOL_HPP
#define DATA_HPOL_HPP

#include<Rcpp.h>
#include <cppad/cppad.hpp>
#include"lgamma.hpp"

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


struct Data {
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

  explicit Data(const Rcpp::List& data)
  : QsansDiag( DgTView(Rcpp::S4(data["QsansDiag"])) )  
    , A(         DgCView(Rcpp::S4(data["ATp"])) )
    , X(         DgCView(Rcpp::S4(data["XTp"])) )
    , CC(        DgCView(Rcpp::S4(data["cc_matrixTp"])) )
    , Qdiag(     data["Qdiag"] )
    , y(         data["y"] )
    , map(       data["map"] )
  {
  // Basic sizes
    Rcpp::S4 QsansDiag2(data["QsansDiag"]);
    Rcpp::NumericVector QsansDiagx=QsansDiag2.slot("x");

    Nq      = static_cast<std::size_t>(QsansDiagx.size());
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


// Thread-safe, R-free copies of your inputs.  switch to DataNoR if threadsafe memory is a concern.
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
  bool dense=0;
  Rcpp::List sparsity;
  Rcpp::List groups;
  Rcpp::List group_sparsity;
  double halfLogDetQ;

  explicit Config(const Rcpp::List& cfg)
  : list(cfg),
  verbose(get_bool(cfg, "verbose", false)),
  dirichlet(get_bool(cfg, "dirichlet", false)),
  transform_theta(get_bool(cfg, "transform_theta", false)),
  beta(get_numvec_copy(cfg, "beta")),
  theta(get_numvec_copy(cfg, "theta")),
  num_threads(get_int(cfg, "num_threads", 1)),           
  dense(get_bool(cfg, "dense", 0)),
  sparsity(cfg.containsElementNamed("sparsity") ? cfg["sparsity"] : Rcpp::List()),
  groups(cfg.containsElementNamed("groups") ? cfg["groups"] : Rcpp::List()),
  group_sparsity(cfg.containsElementNamed("group_sparsity") ? cfg["group_sparsity"] : Rcpp::List()),
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
  Type logNuSq;
  Type oneOverNuSq;
  Type lgammaOneOverNuSq;
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

    int max_val = data.map[0];
    for (std::size_t i = 1; i < data.map.length(); ++i) {
        if (data.map[i] > max_val) {
            max_val = data.map[i];
        }
    }
    Ntheta_base = static_cast<std::size_t>(max_val) + 1;


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
    out.logNuSq           = 2 * out.logTheta[out.logTheta.size()-1];
    out.oneOverNuSq       = CppAD::exp(-out.logNuSq);
    out.lgammaOneOverNuSq = lgamma_any(out.oneOverNuSq);
  } else {
    out.logNuSq = Type(0);
    out.oneOverNuSq = Type(0);
    out.lgammaOneOverNuSq = Type(0);
  }
  return out;
}

inline PackedParams<double>
unpack_params(const std::vector<double>& params,
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

#endif