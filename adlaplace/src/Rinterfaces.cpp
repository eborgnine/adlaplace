#include <Rcpp.h>

Rcpp::List getAdFun_backend(
	Rcpp::List data, 
	Rcpp::List config);
double funH_backend(
	const Rcpp::NumericVector& x,
	const int i,
	SEXP& adPack);
Rcpp::NumericVector gradH_backend(
	const Rcpp::NumericVector& x,
	const int i,
	SEXP  adPack,
	const Rcpp::IntegerVector &pattern);
Rcpp::NumericVector hessH_backend(
  const Rcpp::NumericVector& x,
  const int i,
  SEXP adPack,
  const Rcpp::IntegerVector& row_index,   
  const Rcpp::IntegerVector& col_index  
);

//' C++ backend entry points
//'
//' Low-level C++ entry points exposed to R via Rcpp. These create and operate on
//' an opaque AD “pack” (external pointer) used to evaluate the objective,
//' sparse gradient values, and sparse Hessian values for a selected group.
//'
//' Indices passed in `pattern`, `row`, and `col` are **0-based** and must be
//' sorted with no duplicates (for gradients) or paired consistently (for Hessians).
//'
//' @param data An R list containing model data objects required by the backend.
//' @param config An R list of configuration options required by the backend.
//' @param parameters Numeric vector of parameters.
//' @param i Integer index selecting which group/tape to evaluate.
//' @param adPack External pointer (`externalptr`) returned by \code{getAdFun()}.
//' @param pattern Integer vector of 0-based column indices specifying which
//'   gradient entries to compute (sorted, unique).
//' @param row Integer vector of 0-based row indices for Hessian entries.
//' @param col Integer vector of 0-based column indices for Hessian entries;
//'   must have the same length as \code{row}.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun}: a list containing an opaque external pointer and
//'     associated metadata.
//'   \item \code{jointLogDens}: a scalar objective value for group \code{i}.
//'   \item \code{grad}: numeric vector of length \code{length(pattern)} with
//'     gradient values in the same order as \code{pattern}.
//'   \item \code{hess}: numeric vector of length \code{length(row)} with Hessian
//'     values in the same order as the \code{(row, col)} pairs.
//' }
//'
//' @details
//' The external pointer returned by \code{getAdFun} is opaque and not
//' user-modifiable. It may hold substantial memory (AD tapes, sparsity patterns,
//' work caches). Do not save it across R sessions.
//'
//' @name adlaplace_cpp



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List getAdFun(
	Rcpp::List data, 
	Rcpp::List config)
{

	Rcpp::List result = getAdFun_backend(data, config);

	return result;
}


//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(
	const Rcpp::NumericVector parameters, 
	const int i,
	SEXP adPack)
{
	double result = funH_backend(parameters, i, adPack);

	return(result);
}


//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
	const Rcpp::NumericVector& parameters, 
	const int i,
	SEXP adPack,
	const Rcpp::IntegerVector& pattern)
{

	Rcpp::NumericVector result = gradH_backend(parameters, i, adPack, pattern);

	return(result);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector hess(
	const Rcpp::NumericVector& parameters, 
	const int i,
	SEXP adPack,
	const Rcpp::IntegerVector &row,
	const Rcpp::IntegerVector &col)
{

	Rcpp::NumericVector result = hessH_backend(parameters, i, adPack, row, col);

	return(result);

}
