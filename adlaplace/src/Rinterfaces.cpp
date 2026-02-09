#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/creators/R_interfaces.hpp"

//' C++ backend entry points
//'
//' Low-level C++ entry points exposed to R via Rcpp.
//' These create and operate on an opaque backend handle (external pointer)
//' used to evaluate objective, gradient, and Hessian values.
//'
//' @param data An R list containing model data objects required by the backend
//'   (used by \code{getAdFun()}).
//' @param config An R list of configuration options required by the backend
//'   (used by \code{getAdFun()}).
//' @param x Numeric parameter vector of length \code{Nparams}.
//' @param backendContext External pointer returned by \code{getAdFun()}.
//' @param inner Logical scalar. If \code{TRUE}, evaluate inner-\eqn{\gamma}
//'   derivatives; otherwise evaluate outer/full derivatives.
//' @param Sgroups Optional integer vector of 0-based group indices to evaluate.
//'   If omitted, uses all groups \code{0:(Ngroups-1)}.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun}: external pointer handle with backend state.
//'   \item \code{jointLogDens}: scalar objective value summed over groups.
//'   \item \code{grad}: numeric gradient vector.
//'   \item \code{hess}: sparse symmetric Hessian as a Matrix
//'     \code{dsCMatrix} object.
//' }
//'
//' @details
//' The external pointer returned by \code{getAdFun()} is opaque and not
//' user-modifiable. It may hold substantial memory (AD tapes, sparsity maps,
//' work caches). Do not save it across R sessions.
//'
//' @name adlaplace_cpp



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
SEXP getAdFun(
	Rcpp::List data, 
	Rcpp::List config)
{

	SEXP result = getAdFun_h(data, config);

	return result;
}



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(
	const Rcpp::NumericVector& x,
	SEXP backendContext,
	SEXP Sgroups = R_NilValue) {

	const double fg = jointLogDens_h(x, backendContext, Sgroups);

	return fg;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
	const Rcpp::NumericVector& x,
	SEXP backendContext,
	const bool inner = false,
	SEXP Sgroups = R_NilValue) {

	const Rcpp::NumericVector out = grad_h(x, backendContext, inner, Sgroups);

	return out;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::S4 hess(
	const Rcpp::NumericVector& x,
	SEXP backendContext,
	const bool inner = false,
	SEXP Sgroups = R_NilValue) {

	const Rcpp::S4 out = hess_h(x, backendContext, inner, Sgroups);

	return out;
}



