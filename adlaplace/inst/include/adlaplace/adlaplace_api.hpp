#ifndef ADLAPLACE_API_HPP
#define ADLAPLACE_API_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

double fval_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const size_t i, 
	SEXP adPack);

double grad_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const size_t i, 
	SEXP adPack,
	const bool inner,
	CPPAD_TESTVECTOR(double) &result 
	);

double hess_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const size_t i, 
	SEXP adPack, 
	const bool inner,
	CPPAD_TESTVECTOR(double) &result_grad,
	CPPAD_TESTVECTOR(double) &result_hess
	);

CPPAD_TESTVECTOR(double) thirdDirection_api(
	const CPPAD_TESTVECTOR(double)&  x,
	const CPPAD_TESTVECTOR(double)&  direction,
	const CPPAD_TESTVECTOR(double)&  direction2,  // all zeros
	const CPPAD_TESTVECTOR(double)&  w, // {0.0, 0.0, 1.0}
	const size_t i, 
	SEXP adPack);

#endif

