#include "adlaplace/adlaplace.hpp"

// [[Rcpp::depends(RcppEigen)]]

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


