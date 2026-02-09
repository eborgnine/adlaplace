#include <numeric>   // std::iota

#include "Rcpp.h"
#include "adlaplace/rviews.hpp"


bool get_bool(const Rcpp::List& cfg, const char* key, bool def) {
  return cfg.containsElementNamed(key) ? Rcpp::as<bool>(cfg[key]) : def;
}

int get_int(const Rcpp::List& cfg, const char* key, int def) {
  return cfg.containsElementNamed(key) ? Rcpp::as<int>(cfg[key]) : def;
}

std::vector<int> as_int_vec(
  const Rcpp::IntegerVector& v
  ) {
  std::vector<int> out(static_cast<std::size_t>(v.size()));

  for (R_xlen_t k = 0; k < v.size(); ++k) {
    int x = v[k];
    out[static_cast<std::size_t>(k)] = x;
  }
  return out;
}


CPPAD_TESTVECTOR(double) as_cppad_vector(
  const Rcpp::NumericVector& v
  ) {
  const size_t n = v.size();
  CPPAD_TESTVECTOR(double) out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = v[i];
  }
  return out;
}

// ---- DgCView ----
DgCView::DgCView()
: i(), p(), x(), Dim(Rcpp::IntegerVector::create(0, 0)), has_x(false)
{}

DgCView::DgCView(const Rcpp::S4& obj)
: i(obj.slot("i")),
p(obj.slot("p")),
Dim(obj.slot("Dim")),
has_x(false)
{
  if (obj.inherits("ngCMatrix")) {
    x = Rcpp::NumericVector();   // empty
    has_x = false;
  } else {
    x = Rcpp::as<Rcpp::NumericVector>(obj.slot("x"));
    has_x = true;
  }
}

int DgCView::nrow() const { return Dim[0]; }
int DgCView::ncol() const { return Dim[1]; }
R_xlen_t DgCView::nnz() const { return i.size(); }

// ---- NumVecView ----
NumVecView::NumVecView() : sexp(R_NilValue), is_int(false) {}

NumVecView::NumVecView(SEXP x_) : sexp(x_), is_int(false) {
  const int t = TYPEOF(sexp);
  if (t == INTSXP) {
    is_int = true;
  } else if (t == REALSXP) {
    is_int = false;
  } else {
    Rcpp::stop("NumVecView: expected integer or numeric vector");
  }
}

std::size_t NumVecView::size() const {
  return static_cast<std::size_t>(XLENGTH(sexp));
}

double NumVecView::operator[](std::size_t i_) const {
  return is_int ? static_cast<double>(INTEGER(sexp)[i_]) : REAL(sexp)[i_];
}

// ---- CscPattern ----
CscPattern::CscPattern() : dim(2,0) {}

CscPattern::CscPattern(const Rcpp::S4& sm) : dim(2,0) {
  Rcpp::IntegerVector Ii  = sm.slot("i");
  Rcpp::IntegerVector Pp  = sm.slot("p");
  Rcpp::IntegerVector Dim = sm.slot("Dim");

  dim = as_int_vec(Dim);
  i = as_int_vec(Ii);
  p = as_int_vec(Pp);

  if (sm.inherits("ngCMatrix")) {
    x.resize(i.size());
    std::iota(x.begin(), x.end(), std::size_t(0));
  } else {
    // if this is a dgCMatrix with numeric x, you *can't* read slot("x") as IntegerVector
    // so either:
    //   (A) require x is already integer-coded on the R side, or
    //   (B) accept NumericVector here and cast with checks.
    //
    // You said you're using patterns; so simplest: require ngCMatrix for patterns.
    Rcpp::NumericVector xR = sm.slot("x");
    Rcpp::NumericVector xRound = Rcpp::round(xR, 0);
    Rcpp::IntegerVector xRint = Rcpp::as<Rcpp::IntegerVector>(xRound);
    x = as_int_vec(xRint);
  }
}

/*
MatchGroup::MatchGroup() = default;

MatchGroup::MatchGroup(const Rcpp::List& obj) {
  grad         = CscPattern(Rcpp::as<Rcpp::S4>(obj["grad"]));
  grad_inner   = CscPattern(Rcpp::as<Rcpp::S4>(obj["grad_inner"]));
  hessian      = CscPattern(Rcpp::as<Rcpp::S4>(obj["hessian"]));
  hessian_inner= CscPattern(Rcpp::as<Rcpp::S4>(obj["hessian_inner"]));
}
*/

// ---- Config ----
Config::Config(const Rcpp::List& cfg)
: verbose(get_bool(cfg, "verbose", false)),
transform_theta(get_bool(cfg, "transform_theta", false)),
num_threads(get_int(cfg, "num_threads", 1))
{

  Rcpp::NumericVector beta_nv = cfg["beta"];
  std::vector<double> beta(beta_nv.begin(), beta_nv.end());

  Rcpp::NumericVector gamma_nv = cfg["gamma"];
  std::vector<double> gamma(gamma_nv.begin(), gamma_nv.end());

  Rcpp::NumericVector theta_nv = cfg["theta"];
  std::vector<double> theta(theta_nv.begin(), theta_nv.end());

  beta_begin = 0;
  Nbeta = static_cast<std::size_t>(beta.size());
  beta_end = Nbeta;

  gamma_begin = beta_end;
  Ngamma = static_cast<std::size_t>(gamma.size());
  gamma_end = gamma_begin + Ngamma;

  theta_begin = gamma_end;
  Ntheta = static_cast<std::size_t>(theta.size());
  theta_end = theta_begin + Ntheta;

  Nparams = Nbeta + Ngamma + Ntheta;

  Sgamma.resize(Ngamma);
  params.resize(Nparams);

  for (std::size_t d = 0; d < Nbeta; ++d) {
    params[d] = beta[d];
  }
  for (std::size_t d = 0, idx = gamma_begin; d < Ngamma; ++d, ++idx) {
    Sgamma[d] = (int) idx;
    params[idx] = gamma[d];
  }
  for (std::size_t d = 0, idx = theta_begin; d < Ntheta; ++d, ++idx) {
    params[idx] = theta[d];
  }

  Ngroups = 1;
  if (cfg.containsElementNamed("groups")) {
    groups = CscPattern(Rcpp::as<Rcpp::S4>(cfg["groups"]));
    Ngroups = groups.ncol();
  }

  if (cfg.containsElementNamed("Sgroups")) {
    Sgroups = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(cfg["Sgroups"]));
  } else {
    Sgroups.resize(Ngroups);
    std::iota(Sgroups.begin(), Sgroups.end(), std::size_t(0));
  }
}

std::array<HessianPack,2> hessianPackFromList(const Rcpp::List &x) {

  std::array<HessianPack,2> result;
  HessianPack &hessian_inner = result[0];
  HessianPack &hessian_outer = result[1];


  if (! x.containsElementNamed("hessian")) {
    Rcpp::Rcout << "hessians missing, input should be list(hessian =list(inner=... outer=..))\n";
    return result;
  }

  const Rcpp::List hessian = x["hessian"];

  if (! hessian.containsElementNamed("outer") || ! hessian.containsElementNamed("inner")) {
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
 hessian_inner.map_p = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(inner["p"]));
 hessian_inner.map_local = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(inner["local"]));
 hessian_inner.map_global  = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(inner["global"]));


 Rcpp::List outer = map["outer"];
 hessian_outer.map_p  = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(outer["p"]));
 hessian_outer.map_local   = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(outer["local"]));
 hessian_outer.map_global  = as_int_vec(Rcpp::as<Rcpp::IntegerVector>(outer["global"]));

 return result;
}



Data::Data(const Rcpp::List& data)
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
#ifdef DEBUG    
    if(Ny != A.ncol() && A.nrow() != 0) {
      Rcpp::Rcout << "Ny " << Ny << " columns of A " << A.ncol() << "\n";
    }
    if(Ny != X.ncol() && X.nrow() != 0) {
      Rcpp::Rcout << "lengh y " << Ny << " columns of X " << X.ncol() << "\n";
    }
#endif    
  }