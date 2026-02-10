#ifndef ADLAPLACE_DEFS_HPP
#define ADLAPLACE_DEFS_HPP

#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include <numeric>

#include "adlaplace/runtime/backend.hpp"


// Lightweight view into Matrix::*gCMatrix slots (points to R memory)
struct DgCView {
  Rcpp::IntegerVector i;
  Rcpp::IntegerVector p;
  Rcpp::NumericVector x;     // empty if ngCMatrix
  Rcpp::IntegerVector Dim;
  bool has_x;

  DgCView();
  explicit DgCView(const Rcpp::S4& obj);

  int nrow() const;
  int ncol() const;
  R_xlen_t nnz() const;

  // template must live in the header
  template <class T = double>
  T value(R_xlen_t k) const {
    return has_x ? static_cast<T>(x[k]) : static_cast<T>(1);
  }
};

struct NumVecView {
  SEXP sexp;
  bool is_int;

  NumVecView();
  explicit NumVecView(SEXP x);

  std::size_t size() const;
  double operator[](std::size_t i) const;
};

// Thread-safe CSC pattern copied into std::vector
struct CscPattern {
  std::vector<int> i;
  std::vector<int> p;
  std::vector<int> x;   // if ngCMatrix, filled with 0..nnz-1
  std::vector<int> dim; 
  
  int nrow() const{ return dim[0]; }
  int ncol() const{ return dim[1]; }

  CscPattern();
  explicit CscPattern(const Rcpp::S4& sm);

  std::size_t nnz() const { return i.size(); }
};



struct Config {
  bool verbose;
  bool transform_theta;
  int num_threads;

  std::vector<double> beta, gamma, theta;

  std::size_t beta_begin, beta_end, Nbeta;
  std::size_t gamma_begin, gamma_end, Ngamma;
  std::size_t theta_begin, theta_end, Ntheta;
  std::size_t Nparams, Ngroups;

  std::vector<int> Sgroups;
  std::vector<int> Sgamma;
  std::vector<double> params;

  CscPattern groups;

  explicit Config(const Rcpp::List& cfg);
};


struct Data {
  DgCView A;
  DgCView X;

  NumVecView Qdiag, y;

  DgCView map;
  DgCView elgm_matrix;

  size_t Nmap, Nbeta, Ngamma, Ny;

  explicit Data(const Rcpp::List& data);
};


std::array<HessianPack,2> hessianPackFromList(const Rcpp::List &x);

inline bool adlaplace_get_bool(const Rcpp::List& cfg, const char* key, bool def) {
  return cfg.containsElementNamed(key) ? Rcpp::as<bool>(cfg[key]) : def;
}

inline int adlaplace_get_int(const Rcpp::List& cfg, const char* key, int def) {
  return cfg.containsElementNamed(key) ? Rcpp::as<int>(cfg[key]) : def;
}

inline std::vector<int> adlaplace_as_int_vec(const Rcpp::IntegerVector& v) {
  std::vector<int> out(static_cast<std::size_t>(v.size()));
  for (R_xlen_t k = 0; k < v.size(); ++k) {
    out[static_cast<std::size_t>(k)] = v[k];
  }
  return out;
}

inline DgCView::DgCView()
  : i(), p(), x(), Dim(Rcpp::IntegerVector::create(0, 0)), has_x(false)
{}

inline DgCView::DgCView(const Rcpp::S4& obj)
  : i(obj.slot("i")),
    p(obj.slot("p")),
    Dim(obj.slot("Dim")),
    has_x(false)
{
  if (obj.inherits("ngCMatrix")) {
    x = Rcpp::NumericVector();
    has_x = false;
  } else {
    x = Rcpp::as<Rcpp::NumericVector>(obj.slot("x"));
    has_x = true;
  }
}

inline int DgCView::nrow() const { return Dim[0]; }
inline int DgCView::ncol() const { return Dim[1]; }
inline R_xlen_t DgCView::nnz() const { return i.size(); }

inline NumVecView::NumVecView() : sexp(R_NilValue), is_int(false) {}

inline NumVecView::NumVecView(SEXP x_) : sexp(x_), is_int(false) {
  const int t = TYPEOF(sexp);
  if (t == INTSXP) {
    is_int = true;
  } else if (t == REALSXP) {
    is_int = false;
  } else {
    Rcpp::stop("NumVecView: expected integer or numeric vector");
  }
}

inline std::size_t NumVecView::size() const {
  return static_cast<std::size_t>(XLENGTH(sexp));
}

inline double NumVecView::operator[](std::size_t i_) const {
  return is_int ? static_cast<double>(INTEGER(sexp)[i_]) : REAL(sexp)[i_];
}

inline CscPattern::CscPattern() : dim(2, 0) {}

inline CscPattern::CscPattern(const Rcpp::S4& sm) : dim(2, 0) {
  Rcpp::IntegerVector Ii = sm.slot("i");
  Rcpp::IntegerVector Pp = sm.slot("p");
  Rcpp::IntegerVector Dim = sm.slot("Dim");

  dim = adlaplace_as_int_vec(Dim);
  i = adlaplace_as_int_vec(Ii);
  p = adlaplace_as_int_vec(Pp);

  if (sm.inherits("ngCMatrix")) {
    x.resize(i.size());
    std::iota(x.begin(), x.end(), std::size_t(0));
  } else {
    Rcpp::NumericVector xR = sm.slot("x");
    Rcpp::NumericVector xRound = Rcpp::round(xR, 0);
    Rcpp::IntegerVector xRint = Rcpp::as<Rcpp::IntegerVector>(xRound);
    x = adlaplace_as_int_vec(xRint);
  }
}

inline Config::Config(const Rcpp::List& cfg)
  : verbose(adlaplace_get_bool(cfg, "verbose", false)),
    transform_theta(adlaplace_get_bool(cfg, "transform_theta", false)),
    num_threads(adlaplace_get_int(cfg, "num_threads", 1))
{
  Rcpp::NumericVector beta_nv = cfg["beta"];
  this->beta.assign(beta_nv.begin(), beta_nv.end());

  Rcpp::NumericVector gamma_nv = cfg["gamma"];
  this->gamma.assign(gamma_nv.begin(), gamma_nv.end());

  Rcpp::NumericVector theta_nv = cfg["theta"];
  this->theta.assign(theta_nv.begin(), theta_nv.end());

  beta_begin = 0;
  Nbeta = static_cast<std::size_t>(this->beta.size());
  beta_end = Nbeta;

  gamma_begin = beta_end;
  Ngamma = static_cast<std::size_t>(this->gamma.size());
  gamma_end = gamma_begin + Ngamma;

  theta_begin = gamma_end;
  Ntheta = static_cast<std::size_t>(this->theta.size());
  theta_end = theta_begin + Ntheta;

  Nparams = Nbeta + Ngamma + Ntheta;

  Sgamma.resize(Ngamma);
  params.resize(Nparams);

  for (std::size_t d = 0; d < Nbeta; ++d) {
    params[d] = this->beta[d];
  }
  for (std::size_t d = 0, idx = gamma_begin; d < Ngamma; ++d, ++idx) {
    Sgamma[d] = static_cast<int>(idx);
    params[idx] = this->gamma[d];
  }
  for (std::size_t d = 0, idx = theta_begin; d < Ntheta; ++d, ++idx) {
    params[idx] = this->theta[d];
  }

  Ngroups = 1;
  if (cfg.containsElementNamed("groups")) {
    groups = CscPattern(Rcpp::as<Rcpp::S4>(cfg["groups"]));
    Ngroups = groups.ncol();
  }

  if (cfg.containsElementNamed("Sgroups")) {
    Sgroups = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(cfg["Sgroups"]));
  } else {
    Sgroups.resize(Ngroups);
    std::iota(Sgroups.begin(), Sgroups.end(), std::size_t(0));
  }
}

inline std::array<HessianPack, 2> hessianPackFromList(const Rcpp::List& x) {
  std::array<HessianPack, 2> result;
  HessianPack& hessian_inner = result[0];
  HessianPack& hessian_outer = result[1];

  if (!x.containsElementNamed("hessian")) {
    Rcpp::Rcout << "hessians missing, input should be list(hessian =list(inner=... outer=..))\n";
    return result;
  }

  const Rcpp::List hessian = x["hessian"];
  if (!hessian.containsElementNamed("outer") || !hessian.containsElementNamed("inner")) {
    Rcpp::Rcout << "inner, outer hessian missing, input should be list(hessian =list(inner=... outer=..))\n";
    return result;
  }

  {
    CscPattern hessianCSC(Rcpp::as<Rcpp::S4>(hessian["outer"]));
    hessian_outer.hessian_p = hessianCSC.p;
    hessian_outer.hessian_i = hessianCSC.i;
    hessian_outer.dim = hessianCSC.dim;
  }
  {
    CscPattern hessianCSC(Rcpp::as<Rcpp::S4>(hessian["inner"]));
    hessian_inner.hessian_p = hessianCSC.p;
    hessian_inner.hessian_i = hessianCSC.i;
    hessian_inner.dim = hessianCSC.dim;
  }

  if (!x.containsElementNamed("map")) {
    Rcpp::Rcout << "map missing, input should be list(map =list(inner=... outer=..))\n";
    return result;
  }
  Rcpp::List map = x["map"];
  if (!map.containsElementNamed("inner") || !map.containsElementNamed("outer")) {
    Rcpp::Rcout << "map inner, outer missing, input should be list(map =list(inner=... outer=..))\n";
    return result;
  }

  Rcpp::List inner = map["inner"];
  hessian_inner.map_p = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(inner["p"]));
  hessian_inner.map_local = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(inner["local"]));
  hessian_inner.map_global = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(inner["global"]));

  Rcpp::List outer = map["outer"];
  hessian_outer.map_p = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(outer["p"]));
  hessian_outer.map_local = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(outer["local"]));
  hessian_outer.map_global = adlaplace_as_int_vec(Rcpp::as<Rcpp::IntegerVector>(outer["global"]));

  return result;
}

inline Data::Data(const Rcpp::List& data)
  : A(DgCView(Rcpp::S4(data["ATp"]))),
    X(DgCView(Rcpp::S4(data["XTp"]))),
    Qdiag(data["Qdiag"]),
    y(data["y"]),
    map(DgCView(Rcpp::S4(data["map"])))
{
  if (data.containsElementNamed("elgm_matrix")) {
    elgm_matrix = DgCView(Rcpp::S4(data["elgm_matrix"]));
  }

  Nmap = map.ncol();
  Nbeta = static_cast<std::size_t>(X.nrow());
  Ngamma = static_cast<std::size_t>(A.nrow());
  Ny = static_cast<std::size_t>(y.size());
}

#endif
