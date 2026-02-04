#ifndef ADFUN_API_HPP
#define ADFUN_API_HPP


// needs logDensExtra, logDensRandom, logDensObs to be defined before including

#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <algorithm>
#include <numeric>
#include <vector>

#include "adlaplace/adlaplace_api.hpp"
#include "adlaplace/utils.hpp"

static const std::string JAC_COLOR  = "cppad";
static const std::string HESS_COLOR = "cppad.symmetric";



template <class SizeVec, class ValVec>
inline ValVec& rcv_val(CppAD::sparse_rcv<SizeVec, ValVec>& rcv)
{
  // CppAD sometimes returns val() as const&, even for non-const rcv
	return const_cast<ValVec&>(rcv.val());
}

inline Rcpp::IntegerVector as_int_vec(
	const CPPAD_TESTVECTOR(size_t)& v
	) {
	Rcpp::IntegerVector out(v.size());
	for (size_t i = 0; i < v.size(); ++i) {
		out[i] = static_cast<int>(v[i]);
	}
	return out;
}

inline void adpack_sparsity(
	const CPPAD_TESTVECTOR(double) &x,
  const std::vector<int> &subset,   // sorted, unique, 0-based
  GroupPack& gp,
  const bool verbose=false
  ) {
	const size_t Nparams = x.size();


	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad;
	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad_inner;

	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_upper;

	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner;
	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner_upper;

  // --- seed for for_jac_sparsity (identity n x n) ---
	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> pattern_in;
	pattern_in.resize(Nparams, Nparams, Nparams);
	for (size_t j = 0; j < Nparams; ++j)
		pattern_in.set(j, j, j);

  // --- for_hes_sparsity selectors ---
	CppAD::vectorBool select_domain(Nparams), select_range(1);
	for (size_t j = 0; j < Nparams; ++j) {
		select_domain[j] = true;
	}
	select_range[0] = true;

	const bool transpose     = false;
	const bool dependency    = false;
	const bool internal_bool = false;

	gp.w.resize(1);
	gp.w[0]=1.0;

  // --- full gradient sparsity (m x n, with m = 1 for scalar range) ---
	gp.fun.for_jac_sparsity(
		pattern_in,
		transpose,
		dependency,
		internal_bool,
		grad
		);

  // --- full Hessian sparsity (n x n) ---
	gp.fun.for_hes_sparsity(
		select_domain,
		select_range,
		internal_bool,
		hessian
		);

  // Sanity: scalar output expected (optional, but helpful)
#ifdef DEBUG
	if (grad.nr() != 1 || grad.nc() != Nparams) {
		Rf_error("Expected grad pattern to be 1 x Nparams (scalar output).");
	}
	if (hessian.nr() != Nparams || hessian.nc() != Nparams) {
		Rf_error("Expected hessian pattern to be Nparams x Nparams.");
	}
	if (!std::is_sorted(subset.begin(), subset.end())) {
		Rf_error("subset must be sorted for binary_search.");
	}
#endif


	const auto& full_cols = grad.col();

	std::vector<unsigned char> insubset(Nparams, 0);
	for (size_t v : subset) {
#ifdef DEBUG
		if (v >= Nparams) Rf_error("subset index out of bounds");
#endif
		insubset[v] = 1;             
	}

  // -------- grad_inner = intersection of grad columns with subset --------
	{

// count number in subset
		size_t K = 0;
		for (size_t j : full_cols) {
			K += insubset[j];
		}
		grad_inner.resize(1, Nparams, K);

		size_t t = 0;
		for (size_t j : full_cols) {
			if (insubset[j]) grad_inner.set(t++, 0, j);
		}
	}

  // -------- Hessian derived patterns --------
	{
		const auto& r = hessian.row();
		const auto& c = hessian.col();
		const size_t K = hessian.nnz();

    // count nnz for each derived pattern
		size_t K_upper = 0, K_inner = 0, K_inner_upper = 0;

		for (size_t k = 0; k < K; ++k) {
			const size_t i = r[k], j = c[k];
      const bool is_upper = (i <= j); // choose upper triangle convention
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper) ++K_upper;
      if (is_inner) ++K_inner;
      if (is_upper && is_inner) ++K_inner_upper;
    }

    hessian_upper.resize(Nparams, Nparams, K_upper);
    hessian_inner.resize(Nparams, Nparams, K_inner);
    hessian_inner_upper.resize(Nparams, Nparams, K_inner_upper);

    size_t tu = 0, ti = 0, tiu = 0;
    for (size_t k = 0; k < K; ++k) {
    	const size_t i = r[k], j = c[k];
    	const bool is_upper = (i <= j);
    	const bool is_inner = insubset[i] && insubset[j];

    	if (is_upper) hessian_upper.set(tu++, i, j);
    	if (is_inner) hessian_inner.set(ti++, i, j);
    	if (is_upper && is_inner) hessian_inner_upper.set(tiu++, i, j);
    }
  }

  // evaulate the derivatives in order to set the work 
  gp.fun.Forward(0, x);

  gp.pattern_grad = CppAD::sparse_rcv< CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(grad);
  gp.fun.sparse_jac_rev(
  	x,
  	gp.pattern_grad,
  	grad,
  	JAC_COLOR,      
  	gp.work_grad);

  gp.pattern_grad_inner = CppAD::sparse_rcv< CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double) >(grad_inner);
  gp.fun.sparse_jac_rev(
  	x,
  	gp.pattern_grad_inner,
  	grad_inner,
  	JAC_COLOR,      
  	gp.work_inner_grad);


  gp.pattern_hessian = CppAD::sparse_rcv< CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double) >(hessian_upper);
  gp.fun.sparse_hes(
  	x, gp.w,      
  	gp.pattern_hessian,              
  	hessian,        
  	HESS_COLOR,                
  	gp.work_hess);

  gp.pattern_hessian_inner = CppAD::sparse_rcv< CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double) >(hessian_inner_upper);;
  gp.fun.sparse_hes(
  	x, gp.w,      
  	gp.pattern_hessian_inner,              
  	hessian_inner,        
  	HESS_COLOR,                
  	gp.work_inner_hess);

}


inline std::vector<GroupPack> getAdFun(
	const Data& data,
	const Config& config) {

	size_t NgroupsObs = config.groups.ncol();

	if(NgroupsObs==0) {
		// groups not provided, check for elgm matrix
		NgroupsObs = data.elgm_matrix.ncol();
		if(NgroupsObs == 0) {
			// no elgm matrix, use y
			NgroupsObs = data.y.size();
		}
	}
	const size_t Ngroups = NgroupsObs + 2;

	std::vector<GroupPack> result(Ngroups);	

#ifdef DEBUG
		Rcpp::Rcout << "outer, groups " << Ngroups << " Nbeta " << config.Nbeta << " Ntheta " <<
		config.Ntheta << " Ngamma " << config.Ngamma << " Nparams " << config.Nparams << "\n";
#endif

	CPPAD_TESTVECTOR(double) ad_params_G(config.Nparams);
	for(size_t D=0;D<config.Nparams;D++) {
		ad_params_G[D] = config.params[D];
	}

//# pragma omp parallel
	{
		// local ad_params, in case we're parallel
		CppAD::vector<CppAD::AD<double>> ad_params(config.Nparams);
		for(size_t D=0;D<config.Nparams;D++) {
			ad_params[D] = ad_params_G[D];
		}

//    # pragma omp for schedule(dynamic,1) nowait
		for(size_t D=0;D<NgroupsObs;D++) {

			CppAD::Independent(ad_params);
			auto resultHere = logDensObs(
				ad_params,
				data, config, D);
			CppAD::ADFun<double> fun(ad_params, resultHere);

			result[D].fun = std::move(fun);
			// add sparsity bits
			adpack_sparsity(
				ad_params_G,
				config.Sgamma,
				result[D],
				config.verbose);
		}

//# pragma omp single nowait
		{
			const size_t D = NgroupsObs;
			CppAD::Independent(ad_params);

			auto resultHere = logDensRandom(
				ad_params,
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);
			// add sparsity bits
			adpack_sparsity(
				ad_params_G,
				config.Sgamma,
				result[D],
				config.verbose);
		}

//		# pragma omp single nowait
		{
			const size_t D = NgroupsObs + 1;
			CppAD::Independent(ad_params);

			auto resultHere = logDensExtra(
				ad_params,
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);

			adpack_sparsity(
				ad_params_G,
				config.Sgamma,
				result[D],
				config.verbose);
		}
	} // parallel
	return result;
}

Rcpp::List getAdFun_api(
	const Data& data,
	const Config& config) {

	auto adPack = getAdFun(data, config);

	size_t N = adPack.size();
	Rcpp::List sparsity(N);


	for(size_t D=0;D<N;D++) {
		sparsity[D] = Rcpp::List::create(
			Rcpp::_["grad"] = as_int_vec(adPack[D].pattern_grad.col()), 
			Rcpp::_["grad_inner"] = as_int_vec(adPack[D].pattern_grad_inner.col()),

			Rcpp::_["row_hess"] = as_int_vec(adPack[D].pattern_hessian.row()),
			Rcpp::_["col_hess"] = as_int_vec(adPack[D].pattern_hessian.col()),

			Rcpp::_["row_hess_inner"] = as_int_vec(adPack[D].pattern_hessian_inner.row()),
			Rcpp::_["col_hess_inner"] = as_int_vec(adPack[D].pattern_hessian_inner.col())

			);
	}

	auto* ptr = new std::vector<GroupPack>(std::move(adPack));
	Rcpp::XPtr<std::vector<GroupPack>> xptr(ptr, true);
	xptr.attr("class") = "adpack_ptr";

	Rcpp::List result = Rcpp::List::create(
		Rcpp::_["adPack"] = xptr, 
		Rcpp::_["sparsity"] = sparsity
		);
	return(result);
}

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


CPPAD_TESTVECTOR(double) thirdDirection_api(
	const CPPAD_TESTVECTOR(double)&  x,
	const CPPAD_TESTVECTOR(double)&  direction,
	const CPPAD_TESTVECTOR(double)&  direction2,  // all zeros
	const CPPAD_TESTVECTOR(double)&  w, // {0.0, 0.0, 1.0}
	const size_t i, 
	SEXP adPack
	) {

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];

	gp.fun.Forward(0, x);
	gp.fun.Forward(1, direction);
	gp.fun.Forward(2, direction2);

	return(gp.fun.Reverse(3, w));
}


#endif
