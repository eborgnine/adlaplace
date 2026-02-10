#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/creators/R_interfaces.hpp"
#include "adlaplace/creators/callables.hpp"

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
//'   \item \code{traceHinvT}: numeric vector of third-derivative contractions.
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
  SEXP out = adlaplace_joint_log_dens_callable()(x, backendContext, Sgroups);
  return Rcpp::as<double>(out);
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
  SEXP out = adlaplace_grad_callable()(x, backendContext, Rcpp::wrap(inner), Sgroups);
  return Rcpp::as<Rcpp::NumericVector>(out);
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
  SEXP out = adlaplace_hess_callable()(x, backendContext, Rcpp::wrap(inner), Sgroups);
  return Rcpp::as<Rcpp::S4>(out);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  SEXP backendContext,
  SEXP Sgroups = R_NilValue)
{
  SEXP out = adlaplace_trace_hinv_t_callable()(
    x,
    backendContext,
    LinvPt,
    LinvPtColumns,
    Sgroups
  );
  return Rcpp::as<Rcpp::NumericVector>(out);
}
