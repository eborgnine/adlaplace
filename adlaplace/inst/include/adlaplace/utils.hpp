#ifndef ADLAPLACE_UTILS_HPP
#define ADLAPLACE_UTILS_HPP

#include <Rcpp.h>
#include <vector>
#include <cstddef>


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

struct MatchGroup {
  CscPattern grad, grad_inner, hessian, hessian_inner;

  MatchGroup();
  explicit MatchGroup(const Rcpp::List& obj);
};

struct Config {
  bool verbose;
  bool transform_theta;
  int num_threads;

  Rcpp::NumericVector beta, gamma, theta;

  std::size_t beta_begin, beta_end, Nbeta;
  std::size_t gamma_begin, gamma_end, Ngamma;
  std::size_t theta_begin, theta_end, Ntheta;
  std::size_t Nparams, Ngroups;

  std::vector<int> Sgroups;
  std::vector<int> Sgamma;
  std::vector<double>      params;

  // NOTE: these are patterns; they will copy from S4 into vectors
  CscPattern groups, hessian, hessian_inner, hessian_inner_uplo;

  MatchGroup match;

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

#endif