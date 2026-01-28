
#include "adlaplace/Rinterfaces.hpp"

//' Low-level C++ entry points (Rcpp exports)
//'
//' @param data An R list containing model data and matrices required by the AD
//'   construction. Required fields depend on the model (see vignette).
//' @param config An R list of configuration options (e.g., parameter vectors,
//'   sparsity options, threading options).
//' @param inner Logical; if TRUE, build the inner AD function; otherwise build
//'   the outer AD function.
//'
//' @param parameters Numeric vector of parameters for the requested operation
//'   (e.g., inner optimization or trace calculation). Interpretation depends on
//'   the backend/model.
//' @param control Control list for the optimizer (e.g., trust region settings,
//'   tolerances, max iterations). Used by `inner_opt()`.
//' @param adPackFull,adPack Optional external pointer (`externalptr`) returned by
//'   `getAdFun()`. If provided, reuse cached AD objects / sparsity structures.
//'
//' @param LinvPt,LinvPtColumns Sparse matrix objects (S4, typically from Matrix)
//'   used by `traceHinvT()`; see package documentation/vignette for required
//'   classes and dimensions.
//'
//' @return `getAdFun()` returns an external pointer (`externalptr`) to an
//'   internal AD object. 
//'
//'
//' @details
//' The returned pointer is not human-readable and should not be modified.
//' It may hold substantial memory (tapes, sparsity patterns, caches). Use
//' package-level functions to manage lifecycle and computation.
//'
//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
SEXP getAdFun(
	Rcpp::List data, 
	Rcpp::List config,
	const bool inner=false)
{
	auto xp = getAdFun_backend(data, config, inner);
	return xp;
}



//' @rdname adlaplace_cpp
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

//' @rdname adlaplace_cpp
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



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
	const Rcpp::NumericVector parameters,
	const Rcpp::S4& LinvPt,
	const Rcpp::S4& LinvPtColumns,
	const Rcpp::List config,
	SEXP adPack = R_NilValue
	) {

	auto result = traceHinvT_backend(parameters, 
		LinvPt, LinvPtColumns, config, adPack);
	return(result);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List sparsity(
   Rcpp::List data,
   Rcpp::List config
	) {
	
	auto result=sparsity_backend(data, config);

	return(result);
}


