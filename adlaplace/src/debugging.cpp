

#include "adlaplace/debugging.hpp"

//' Joint log-density and derivatives
//'
//' Compute the joint log-density of a hierarchical model and its first- and
//' second-order derivatives using automatic differentiation.  
//'
//' These functions are low-level computational backends intended for deugging.
//' They operate on a vector of parameters and
//' a pre-constructed automatic differentiation object (`adPack`).
//'
//' @param parameters Numeric vector of model parameters. The ordering must
//'   match the structure expected by the automatic differentiation object.
//' @param adPack An external pointer (SEXP) to an automatic differentiation
//'   object created by \code{getAdFun()}.
//' @param config A list of configuration options controlling model structure
//'   and computation.
//'
//' @return
//' \describe{
//'   \item{\code{jointLogDens()}}{A scalar numeric value giving the joint
//'     log-density.}
//'   \item{\code{jointLogDensNoAdfun()}}{As above, but takes data rather than adPack.}
//'   \item{\code{grad()}}{A numeric vector giving the gradient of the joint
//'     log-density with respect to \code{parameters}.}
//'   \item{\code{hessian()}}{A sparse symmetric matrix (as a
//'     \code{dgCMatrix}) giving the Hessian of the joint log-density.}
//' }
//'
//' @details
//' These functions rely on CppAD-based automatic differentiation. The
//' \code{adPack} object encapsulates the taped computation graph and must be
//' compatible with the supplied \code{parameters} and \code{config}.
//'
//' Users should generally call higher-level R wrappers rather than invoking
//' these functions directly.
//'
//' @seealso \code{\link{getAdFun}}, \code{\link{logLik}}
//'
//' @name jointLogDens
//' @rdname jointLogDens
//' @export
// [[Rcpp::export]]
double jointLogDens(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	double result = jointLogDens_backend(parameters, adPack, config);

	return(result);

}

//' @param data List of data objects used to construct the AD function.
//' @rdname jointLogDens
//' @export
// [[Rcpp::export]]
double jointLogDensNoAdfun(
	Rcpp::NumericVector parameters, 
	Rcpp::List data,
	Rcpp::List config)
{

	double result = jointLogDensNoAdfun_backend(parameters, data, config);
	return(result);
}



//' @rdname jointLogDens
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	auto result = grad_backend(parameters, adPack, config);
	return(result);

}

//' @rdname jointLogDens
//' @export
// [[Rcpp::export]]
Rcpp::S4 hessian(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	auto result = hessian_backend(parameters, adPack, config);

	return(result);
}