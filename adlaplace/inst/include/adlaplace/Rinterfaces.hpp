#ifndef RINTERFACES_HPP
#define RINTERFACES_HPP

// a true header file, includes forward declarations only

#include <Rcpp.h>

// ---- Forward declarations only (no CppAD types here) ----
SEXP getAdFun_backend(Rcpp::List data, Rcpp::List config, bool inner);

Rcpp::List inner_opt_backend(
  Rcpp::NumericVector parameters,
  Rcpp::List data,
  Rcpp::List control,
  Rcpp::List config,
  SEXP adPackFull
);

Rcpp::List inner_opt_adpack_backend(
  Rcpp::NumericVector parameters,
  SEXP adPack,
  Rcpp::List control,
  Rcpp::List config,
  SEXP adPackFull
);

Rcpp::NumericVector traceHinvT_backend(
  Rcpp::NumericVector parameters,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  Rcpp::List config,
  SEXP adPack
);

Rcpp::List sparsity_backend(Rcpp::List data, Rcpp::List config);


double jointLogDens_backend(
	Rcpp::NumericVector, 
	SEXP,
	Rcpp::List);

Rcpp::NumericVector grad_backend(
	Rcpp::NumericVector, 
	SEXP,
	Rcpp::List);

Rcpp::S4 hessian_backend(
	Rcpp::NumericVector, 
	SEXP,
	Rcpp::List);

double jointLogDensNoAdfun_backend(
	Rcpp::NumericVector parameters, 
	Rcpp::List dataR,
	Rcpp::List configR);



#endif

