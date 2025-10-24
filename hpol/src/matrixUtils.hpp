#include "hpol.hpp"

Rcpp::S4 make_TMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N);

Rcpp::S4 make_CMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& p);
Rcpp::S4 make_gCMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& p);

Rcpp::IntegerVector  compute_p_vector(
  const Rcpp::IntegerVector & j, int ncol);

CPPAD_TESTVECTOR( std::set<size_t> ) build_pattern_from_R(
  const Rcpp::IntegerVector& row0,
  const Rcpp::IntegerVector& col0,
  size_t n);

Rcpp::S4 assembleHessian(
      const std::vector<std::vector<double>>& randomHessian, 
      const std::vector<double>& qHessian, 
      const Rcpp::List& sparsity, 
      const Config& config, 
      const bool onlyRandom);
