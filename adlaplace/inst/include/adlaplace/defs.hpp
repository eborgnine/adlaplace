#ifndef ADLAPLACE_DEFS_HPP
#define ADLAPLACE_DEFS_HPP

#include <Rcpp.h>
#include <vector>
#include <cstddef>


struct HessianPack {
  std::vector<int> hessian_p;
  std::vector<int> hessian_i;
  std::vector<int> dim;

  std::vector<int> map_p;
  std::vector<int> map_local;
  std::vector<int> map_global;

};

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

  Rcpp::NumericVector beta, gamma, theta;

  std::size_t beta_begin, beta_end, Nbeta;
  std::size_t gamma_begin, gamma_end, Ngamma;
  std::size_t theta_begin, theta_end, Ntheta;
  std::size_t Nparams, Ngroups;

  std::vector<int> Sgroups;
  std::vector<int> Sgamma;
  std::vector<double> params;

  // NOTE: these are patterns; they will copy from S4 into vectors
  CscPattern groups,  hessian_inner_uplo; 
  HessianPack hessian, hessian_inner;

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




struct GroupPack {
  CppAD::ADFun<double>              fun;       // taped function for the group
  CppAD::sparse_jac_work            work_grad;
  CppAD::sparse_hes_work            work_hess;      // reusable work cache
  CppAD::sparse_jac_work            work_inner_grad;
  CppAD::sparse_hes_work            work_inner_hess;      // reusable work cache
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad_inner;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian;  // note these are upper triangle only
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian_inner;

  CPPAD_TESTVECTOR(double) w;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

  CPPAD_TESTVECTOR(double) x;

};


// backend private context
struct BackendContext {
  std::vector<GroupPack> *adPack;
//  std::array<std::array<int,3>,3> sizes; // beta, gamma, theta; begin end N

  size_t Nparams;
  size_t Ngroups;
  size_t Nbeta;
  size_t Ngamma;
  size_t Ntheta;

  HessianPack hessian_inner; // p, local index, global index, dims
  HessianPack hessian_outer; // 

};


#endif