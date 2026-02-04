#ifndef ADLAPLACE_API_HPP
#define ADLAPLACE_API_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>


struct GroupPack {
  CppAD::ADFun<double>              fun;       // taped function for the group
  CppAD::sparse_jac_work            work_grad;
  CppAD::sparse_hes_work            work_hess;      // reusable work cache
  CppAD::sparse_jac_work            work_inner_grad;
  CppAD::sparse_hes_work            work_inner_hess;      // reusable work cache
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad_inner;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian;  // note these are upper triangle only
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian_inner;

  CPPAD_TESTVECTOR(double) w;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

  CPPAD_TESTVECTOR(double) x;

};

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

