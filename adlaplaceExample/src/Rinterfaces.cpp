#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/creators/R_interfaces.hpp"

//' C++ backend entry points for adlaplaceExample
//'
//' Build and evaluate skew-normal model AD handles compatible with
//' \pkg{adlaplace} runtime routines.
//'
//' @param data Model data list.
//' @param config Model configuration list.
//' @param x Full parameter vector.
//' @param backendContext External pointer or list containing \code{adFun}.
//' @param inner Logical; if \code{TRUE}, evaluate inner-parameter derivatives.
//' @param Sgroups Optional 0-based integer group indices. Defaults to all groups.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun}: list with \code{adFun}, \code{sparsity}, and \code{hessians}.
//'   \item \code{jointLogDens}: scalar objective value.
//'   \item \code{grad}: numeric gradient vector.
//'   \item \code{hess}: sparse symmetric Hessian as \code{dsCMatrix}.
//' }
//'
//' @name adlaplace_cpp

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List getAdFun(
  Rcpp::List data,
  Rcpp::List config)
{
  return getAdFun_h(data, config);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(
  const Rcpp::NumericVector& x,
  SEXP backendContext,
  SEXP Sgroups = R_NilValue)
{
  return jointLogDens_h(x, backendContext, Sgroups);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
  const Rcpp::NumericVector& x,
  SEXP backendContext,
  const bool inner = false,
  SEXP Sgroups = R_NilValue)
{
  return grad_h(x, backendContext, inner, Sgroups);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::S4 hess(
  const Rcpp::NumericVector& x,
  SEXP backendContext,
  const bool inner = false,
  SEXP Sgroups = R_NilValue)
{
  return hess_h(x, backendContext, inner, Sgroups);
}
