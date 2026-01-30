#ifndef ADFUN_API_HPP
#define ADFUN_API_HPP

// needs logDensExtra, logDensRandom, logDensObs

// needs Rcpp.hpp, cppad.hpp,  

#include <Rcpp.h>

#include "adlaplace/data.hpp"
#include "adlaplace/foromp.hpp"
#include "adlaplace/cppadUtils.hpp"

inline void adpack_sparsity(
	const CPPAD_TESTVECTOR(double) &x,
	PatternPair& pattern,
	GroupPack& gp,
	const bool verbose=false
	) {

	const size_t Nparams = x.size();
	CPPAD_TESTVECTOR(double) w(1);
	w[0] = 1.0;

	CppAD::vectorBool select_domain(Nparams), select_range(1);
	CppAD::sparse_rc< CPPAD_TESTVECTOR(size_t) > pattern_in;
    pattern_in.resize(Nparams, Nparams, Nparams);   // n nonzeros
    for (size_t j = 0; j < Nparams; ++j) {
    	select_domain[j] = true;
      	pattern_in.set(j, j, j);   // (row=j, col=j, index=j)
      }

    select_range[0] = true; // scalar output
    bool transpose    = false;   // want pattern_out as m × n
    bool dependency   = false;   // standard Jacobian sparsity (not dependency)
    bool internal_bool = false;  // let CppAD choose representation; updated on return

    gp.fun.for_jac_sparsity(
    	pattern_in,
    	transpose,
    	dependency,
    	internal_bool,
    	pattern.first);

    gp.fun.for_hes_sparsity(
    	select_domain,
    	select_range,
    	internal_bool,
    	pattern.second);


    CppAD::sparse_rcv< CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double) > out_grad(pattern.first);
    gp.fun.sparse_jac_rev(
    	x,
    	out_grad,
    	pattern.first,
    	JAC_COLOR,      
    	gp.work_grad);

    CppAD::sparse_rcv< CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double) > out_hess(pattern.second);
    gp.fun.sparse_hes(
    	x, w,      
    	out_hess,              
    	pattern.second,        
    	HESS_COLOR,                
    	gp.work_hess);
}


inline void getAdFun(
	std::vector<GroupPack> &result,
	PatternPairs & pattern,
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

	result.clear();
	result.resize(Ngroups);
	pattern.clear();
	pattern.resize(Ngroups);
	

	const size_t Ngamma = config.start_gamma.size();
	const size_t Nbeta = config.beta.size();
	const size_t Ntheta = config.theta.size();
	const size_t Nparams = Nbeta + Ngamma + Ntheta;

	if(config.verbose) {
		Rcpp::Rcout << "outer, groups " << Ngroups << " Nbeta " << Nbeta << " Ntheta " <<
		Ntheta << " Ngamma " << Ngamma << " Nparams " << Nparams << "\n";
	}

	CPPAD_TESTVECTOR(double) ad_params_G(Nparams);
	for(size_t D=0;D<Nbeta;D++) {
		ad_params_G[D] = config.beta[D];
	}
	for(size_t Dgamma=0;Dgamma<Ngamma;Dgamma++) {
		ad_params_G[Dgamma+Nbeta] = config.start_gamma[Dgamma];
	}
	for(size_t Dtheta=0;Dtheta<Ntheta;Dtheta++) {
		ad_params_G[Dtheta+Nbeta+Ngamma] = config.theta[Dtheta];
	} 

//# pragma omp parallel
	{
		CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
		for(size_t D=0;D<Nparams;D++) {
			ad_params[D] = ad_params_G[D];
		}

//    # pragma omp for schedule(dynamic,1) nowait
		for(size_t D=0;D<NgroupsObs;D++) {


			CppAD::Independent(ad_params);
			auto resultHere = logDensObs(
				slice(ad_params, Nbeta, Nbeta + Ngamma), // gamma
				slice(ad_params, 0, Nbeta), // beta
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config, D);
			CppAD::ADFun<double> fun(ad_params, resultHere);

			result[D].fun = std::move(fun);
			// add sparsity bits
			adpack_sparsity(
				ad_params_G,
				pattern[D],
				result[D],
				config.verbose);
		}

//# pragma omp single nowait
		{
			CppAD::Independent(ad_params);

			auto resultHere = logDensRandom(
				slice(ad_params, Nbeta, Nbeta + Ngamma), // gamma
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[NgroupsObs].fun = std::move(fun);
			// add sparsity bits
			adpack_sparsity(
				ad_params_G,
				pattern[NgroupsObs],
				result[NgroupsObs],
				config.verbose);
		}

//		# pragma omp single nowait
		{
			CppAD::Independent(ad_params);

			auto resultHere = logDensExtra(
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[NgroupsObs+1].fun = std::move(fun);

			// add sparsity bits
			adpack_sparsity(
				ad_params_G,
				pattern[NgroupsObs+1],
				result[NgroupsObs+1],
				config.verbose);
		}
	} // parallel

}



Rcpp::List getAdFun_backend(
	Rcpp::List data, 
	Rcpp::List config)
{

	Data dataC(data);
	Config configC(config);

//	cppad_parallel_setup(configC.num_threads);

	std::vector<GroupPack> adPack;
	PatternPairs pattern;
	getAdFun(adPack, pattern, dataC, configC);

	auto* ptr = new std::vector<GroupPack>(std::move(adPack));
  	Rcpp::XPtr<std::vector<GroupPack>> xp(ptr, /*deleteOnFinalizer=*/true);
	xp.attr("class") = "adpack_ptr";

	Rcpp::List result = Rcpp::List::create(
		Rcpp::_["adPack"] = xp, 
		Rcpp::_["pattern"] = sparsityList(pattern)
		);

	return result;
}

void funH(
	const CPPAD_TESTVECTOR(double) &x, 
	const int i, 
	SEXP adPack, 
	double& fout) {

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];

	fout = gp.fun.Forward(0, x)[0];
}

double funH_backend(
	const Rcpp::NumericVector& x,
	const int i,
	SEXP& adPack) {

	double fout;
	const size_t Nparams = x.size();
	CPPAD_TESTVECTOR(double) xC(Nparams);
	for(size_t D=0;D<Nparams;++D) {
		xC[D] = x[D];
	}
	funH(xC, i, adPack, fout);
	return(fout);
}

void gradH(
	const CPPAD_TESTVECTOR(double) &x, 
	const int i, 
	SEXP adPack, 
	CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> &gout) {

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];
	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

#ifdef DEBUG
  if (gp.work_grad.empty()) {
  	Rf_error("grad work doesnt contain sparsity pattern");
  }
#endif

	gp.fun.sparse_jac_rev(
		x,
		gout,
		unused_pattern,
		JAC_COLOR,      
		gp.work_grad);
}

Rcpp::NumericVector gradH_backend(
	const Rcpp::NumericVector& x,
	const int i,
	SEXP  adPack,
	const Rcpp::IntegerVector &pattern) {

	const size_t Nparams = x.size();
	CPPAD_TESTVECTOR(double) xC(Nparams);
	for(size_t D=0;D<Nparams;++D) {
		xC[D] = x[D];
	}
	auto patternC = build_pattern_from_R(Nparams, pattern);
	CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> resultC(patternC);

	gradH(xC, i, adPack, resultC);
	 const auto& resultValues = resultC.val();
	const size_t Nresult = resultC.nnz();
	Rcpp::NumericVector result(Nresult);
	for(size_t D=0;D<Nresult;++D) {
		result[D] = resultValues[D];
	}
	return(result);
}

void hessH(
	const CPPAD_TESTVECTOR(double) &x, 
	const int i, 
	SEXP& adPack, 
	CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> &hout) {

	CPPAD_TESTVECTOR(double) w(1);
	w[0] = 1.0;

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];
	CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

	gp.fun.sparse_hes(
		x,  
		w,      
		hout,              
		unused_pattern, // not used        
		HESS_COLOR,                
		gp.work_hess               
		);
}

Rcpp::NumericVector hessH_backend(
  const Rcpp::NumericVector& x,
  const int i,
  SEXP adPack,
  const Rcpp::IntegerVector& row_index,   // row indices (0-based)
  const Rcpp::IntegerVector& col_index   // col indices (0-based)
){
  const size_t Nparams = (size_t) x.size();

  // x -> CPPAD vector
  CPPAD_TESTVECTOR(double) xC(Nparams);
  for (size_t k = 0; k < Nparams; ++k)
    xC[k] = x[(int)k];

  // Build Hessian sparsity pattern (N x N) from (row, col) triplets
  auto patternC = build_pattern_from_R(Nparams, row_index, col_index);

  // Output container aligned with pattern
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> hout(patternC);

  // Evaluate sparse Hessian values
  hessH(xC, i, adPack, hout);

  // Copy nnz values back to R
  const auto& vals = hout.val();
  const size_t nnz = hout.nnz();

  Rcpp::NumericVector result(nnz);

  for (size_t k = 0; k < nnz; ++k)
    result[(int)k] = vals[k];

  return(result);
}

inline void thirdDirection(
	CPPAD_TESTVECTOR(double) & result,
	const CPPAD_TESTVECTOR(double)&  parameters,
	const CPPAD_TESTVECTOR(double)&  direction,
	const CPPAD_TESTVECTOR(double)&  direction2,  // all zeros
	const CPPAD_TESTVECTOR(double)&  w, // {0.0, 0.0, 1.0}
	const int i, 
	SEXP adPack
	) {

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	GroupPack &gp = (*xp)[i];

	gp.fun.Forward(0, parameters);
	gp.fun.Forward(1, direction);
	gp.fun.Forward(2, direction2);

	result = gp.fun.Reverse(3, w);

}



#endif
