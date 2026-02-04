#pragma once

#ifndef ADLAPLACE_HPP
#define ADLAPLACE_HPP

// We do NOT support ColPack in this package.#ifndef CPPAD_HAS_COLPACK
#define CPPAD_HAS_COLPACK 0

#include <cppad/cppad.hpp>

// needs cppad but not Rcpp
#include "adlaplace/cppadUtils.hpp"
#include "adlaplace/logexp.hpp"

#include "adlaplace/adpack.hpp"
#include "adlaplace/foromp.hpp"

// needs Eigen, adpack
#include <Eigen/SparseCore>   
#include "adlaplace/functions.hpp"


// for what follows need rcpp

#include <Rcpp.h>
#include "adlaplace/trustOptimUtils.hpp"
// need Eigen 
#include "adlaplace/matrixUtils.hpp"
// needs matrixUtils (DGCVIEW)
#include "adlaplace/data.hpp"


// needs adpack
#include "adlaplace/sparsity.hpp"

// needs adpack, omp, config
#include "adlaplace/third.hpp"

// needs config, data, sparsity
#include "adlaplace/adfun.hpp"

// needs adfun
#include "adlaplace/debugging.hpp"


// from trustOptim
#include <CG-sparse.h> 

// needs functions.hpp and trustoptim
#include "adlaplace/innerOpt.hpp"

// needs everything
#include "adlaplace/Rinterfaces_backend.hpp"



double fval_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const size_t i, 
	SEXP adPack) {

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];
	
	const double result = gp.fun.Forward(0, x)[0];
	return result;
}

double grad_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const size_t i, 
	SEXP adPack,
	const bool inner,
	CPPAD_TESTVECTOR(double) &result 
	){

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];

	auto* pattern = inner?&gp.pattern_grad_inner:&gp.pattern_grad;
	auto* work = inner ? &gp.work_inner_grad : &gp.work_grad;


#ifdef DEBUG
	if (!inner && pattern->nc() != x.size()) {
		Rcpp::Rcout << "grad inner " << inner << " patternnc " << pattern->nc() << " xsize " <<
		x.size() << "\n";
		Rf_error("grad pattern and parameters different lengths");
	}
	if (pattern->nnz() != result.size()) {
		Rf_error("grad pattern and result different lengths");
	}
#endif

	const double result_f = gp.fun.Forward(0, x)[0];

	rcv_val(*pattern).swap(result);

	gp.fun.sparse_jac_rev(
		x,
		*pattern,
		gp.unused_pattern,
		JAC_COLOR,      
		*work);

	rcv_val(*pattern).swap(result);
	return result_f;
}

double hess_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const size_t i, 
	SEXP adPack, 
	const bool inner,
	CPPAD_TESTVECTOR(double) &result_grad,
	CPPAD_TESTVECTOR(double) &result_hess
	){


	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];

	auto* pattern_hess = inner ? &gp.pattern_hessian_inner : &gp.pattern_hessian;
	auto* work_hess = inner ? &gp.work_inner_hess : &gp.work_hess;

	auto* pattern_grad = inner?&gp.pattern_grad_inner:&gp.pattern_grad;
	auto* work_grad = inner ? &gp.work_inner_grad : &gp.work_grad;

	const double result_f = gp.fun.Forward(0, x)[0];

	if(result_grad.size() > 0) {
		rcv_val(*pattern_grad).swap(result_grad);

		gp.fun.sparse_jac_rev(
			x,
			*pattern_grad,
			gp.unused_pattern,
			JAC_COLOR,
			*work_grad);

		rcv_val(*pattern_grad).swap(result_grad);
	}

	rcv_val(*pattern_hess).swap(result_hess);

	gp.fun.sparse_hes(
		x,  
		gp.w,
		*pattern_hess,              
		gp.unused_pattern, // not used        
		HESS_COLOR,                
		*work_hess              
		);

	rcv_val(*pattern_hess).swap(result_hess);

	return result_f;
}

#endif
