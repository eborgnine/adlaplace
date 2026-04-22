#ifndef ADFUN_CREATE_HPP
#define ADFUN_CREATE_HPP



// to be included in objectiveFunction.cpp to create AD functions.

// needs logDensExtra, logDensRandom, logDensObs to be defined before including


#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <algorithm>
#include <numeric>
#include <vector>

#include "adlaplace/adlaplace.hpp" // data and config


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

Rcpp::List extract_sparsity(
	const std::vector<GroupPack> &adFun) {

	const size_t N = adFun.size();
	Rcpp::List sparsity(N);

	for(size_t D=0;D<N;D++) {
		sparsity[D] = Rcpp::List::create(
			Rcpp::_["grad"] = as_int_vec(adFun[D].pattern_grad.col()), 
			Rcpp::_["grad_inner"] = as_int_vec(adFun[D].pattern_grad_inner.col()),

			Rcpp::_["row_hess"] = as_int_vec(adFun[D].pattern_hessian.row()),
			Rcpp::_["col_hess"] = as_int_vec(adFun[D].pattern_hessian.col()),

			Rcpp::_["row_hess_inner"] = as_int_vec(adFun[D].pattern_hessian_inner.row()),
			Rcpp::_["col_hess_inner"] = as_int_vec(adFun[D].pattern_hessian_inner.col())

			);
	}

	return(sparsity);
}

inline void adpack_sparsity(
	const CPPAD_TESTVECTOR(double) &x,
  const std::vector<int> &subset,   // sorted, unique, 0-based
  GroupPack& gp,
  const bool verbose=false
  ) {
	const size_t Nparams = x.size();

	gp.w.resize(1);
	gp.w[0]=1.0;

	gp.wthree.resize(3);
	gp.wthree[0] = 0.0;
	gp.wthree[1] = 0.0;
	gp.wthree[2] = 1.0;

	gp.x.resize(Nparams);
	gp.direction.resize(Nparams);
	gp.direction_zeros.resize(Nparams);

	std::fill(gp.direction.begin(), gp.direction.end(), 0.0);
	std::fill(gp.direction_zeros.begin(), gp.direction_zeros.end(), 0.0);

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


std::vector<GroupPack> getAdFun(
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

if(config.verbose) {
	Rcpp::Rcout << "outer, groups " << Ngroups << " Nbeta " << config.Nbeta << " Ntheta " <<
	config.Ntheta << " Ngamma " << config.Ngamma << " Nparams " << config.Nparams << "\n";
}

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
if(config.verbose) {
	Rcpp::Rcout << ".";
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
if(config.verbose) {
	Rcpp::Rcout << ".";
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
if(config.verbose) {
	Rcpp::Rcout << ".";
}
//		# pragma omp single nowait
		{
			const size_t D = NgroupsObs + 1;
			CppAD::Independent(ad_params);

if(config.verbose) {
	Rcpp::Rcout << "extra ";
}

			auto resultHere = logDensExtra(
				ad_params,
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);

if(config.verbose) {
	Rcpp::Rcout << ".";
}
			adpack_sparsity(
				ad_params_G,
				config.Sgamma,
				result[D],
				config.verbose);
if(config.verbose) {
	Rcpp::Rcout << ".";
}

		}
	} // parallel
if(config.verbose) {
	Rcpp::Rcout << " done.\n";
}
	return result;
}



#endif
