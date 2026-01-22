#include "adlaplace/adlaplace.hpp"

// [[Rcpp::depends(RcppEigen)]]

/*
 *' Construct and return an AD function object (external pointer)
 *'
 *' Creates (or retrieves) an internal automatic-differentiation (AD) function
 *' object used by \pkg{adlaplace} for evaluating the objective, gradients, and
 *' (optionally sparse) derivative structures. The returned object is an
 *' external pointer that is intended to be passed back to other \pkg{adlaplace}
 *' routines (e.g., inner/outer optimization and derivative evaluation).
 *'
 *' This is a low-level interface. Most users should call higher-level R
 *' functions (e.g., \code{logLik()} or \code{inner_opt()}) rather than using
 *' \code{getAdFun()} directly.
 *'
 *' @param data An R list containing model data and matrices required by the AD
 *'   construction. The required fields depend on your model setup (see package
 *'   vignette).
 *' @param config An R list of configuration options (e.g., parameter vectors,
 *'   sparsity options, threading options). The expected entries depend on the
 *'   model and backend.
 *' @param inner Logical; if \code{TRUE}, build the "inner" AD function (typically
 *'   corresponding to the random-effect subproblem). If \code{FALSE} (default),
 *'   build the "outer" AD function.
 *'
 *' @return An external pointer (\code{externalptr}) to an internal AD object.
 *'   This pointer is meant to be used only by \pkg{adlaplace} functions and is
 *'   not stable across sessions.
 *'
 *' @details
 *' The returned pointer is not human-readable and should not be modified.
 *' It may hold substantial memory (tapes, sparsity patterns, caches). Use
 *' package-level functions to manage lifecycle and computation.
 *'
 *'
 *' @export
 */
// [[Rcpp::export]]
SEXP getAdFun(
	Rcpp::List data, 
	Rcpp::List config,
	const bool inner=false)
{

	auto xp = getAdFun_backend(data, config, inner);
	return xp;
}

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

//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
	Rcpp::NumericVector parameters, 
	Rcpp::List data,
	Rcpp::List control, 
	Rcpp::List config,
	SEXP adPackFull = R_NilValue)
{
	auto res = inner_opt_backend(parameters, data, control, config, adPackFull);
	return(res);
}


//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt_adpack(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	const Rcpp::List control, 
	const Rcpp::List config,
	SEXP adPackFull = R_NilValue)
{


	auto res = inner_opt_adpack_backend(parameters, adPack, control, config, adPackFull);
	return(res);
}



//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
	const Rcpp::NumericVector parameters,
	const Rcpp::S4& LinvPt,
	const Rcpp::S4& LinvPtColumns,
	const Rcpp::List config,
	SEXP adPack = R_NilValue
	) {

	auto result = traceHinvT_backend(parameters, LinvPt, LinvPtColumns, config, adPack);
	return(result);
}

//' @export
// [[Rcpp::export]]
Rcpp::List sparsity(
   Rcpp::List data,
   Rcpp::List config
	) {
	
	auto result=sparsity_backend(data, config);

	return(result);
}


