#ifndef DEBUGGING_HPP
#define DEBUGGING_HPP

#include "adlaplace/adlaplace_base.hpp"

#include "adlaplace/adpack.hpp"
#include "adlaplace/sparsity.hpp"
#include "adlaplace/functions.hpp"
#include "adlaplace/adfun.hpp"

double jointLogDens_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	Config configC(config);
	double result = 0.0;

	if (adPack == R_NilValue) {
		return(result);
	}
    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}
	const bool inner = (Nparams == configC.start_gamma.size());
	auto hessianIPLowerHere = inner?configC.hessianIPLower_inner:configC.hessianIPLower_outer;


	AD_Func_Opt funObj(*xp, hessianIPLowerHere);

	Eigen::VectorXd parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	set_num_threads_wrapper(configC.num_threads);
	CppAD::thread_alloc::parallel_setup(
		configC.num_threads,
		[](){ return in_parallel_wrapper(); },
		[](){ return static_cast<size_t>(thread_num_wrapper()); }
		);

	if(configC.verbose) {
		Rcpp::Rcout << "d";
	}


	funObj.get_f(parametersC, result);

	if(configC.verbose) {
		Rcpp::Rcout << "e " << result << "\n";
	}


	return(result);

}


Rcpp::NumericVector grad_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{
	Config configC(config);

	if (adPack == R_NilValue) {
		return(Rcpp::NumericVector());
	}

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}
	const bool inner = (Nparams == configC.start_gamma.size());
	auto hessianIPLowerHere = inner?configC.hessianIPLower_inner:configC.hessianIPLower_outer;

	AD_Func_Opt funObj(*xp, hessianIPLowerHere);

	Rcpp::NumericVector result(Nparams);


	Eigen::VectorXd parametersC(Nparams);
	Eigen::VectorXd resultC(Nparams);

	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	set_num_threads_wrapper(configC.num_threads);
	CppAD::thread_alloc::parallel_setup(
		configC.num_threads,
		[](){ return in_parallel_wrapper(); },
		[](){ return static_cast<size_t>(thread_num_wrapper()); }
		);

	funObj.get_df(parametersC, resultC);

	for(size_t D=0; D<Nparams;D++) {
		result[D] = resultC[D];
	}
	return(result);

}

Rcpp::S4 hessian_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	if (adPack == R_NilValue) {
		return(R_NilValue);
	}

	Config configC(config);

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}

	const bool inner = (Nparams == configC.start_gamma.size());
	auto hessianIPLowerHere = inner?configC.hessianIPLower_inner:configC.hessianIPLower_outer;

	AD_Func_Opt funObj(*xp, hessianIPLowerHere);


	Eigen::VectorXd parametersC(Nparams);
	
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	Eigen::SparseMatrix<double> resultC = funObj.Htemplate.cast<double>();

	set_num_threads_wrapper(configC.num_threads);
	CppAD::thread_alloc::parallel_setup(
		configC.num_threads,
		[](){ return in_parallel_wrapper(); },
		[](){ return static_cast<size_t>(thread_num_wrapper()); }
		);

	if(configC.verbose) {
		Rcpp::Rcout << "hessian params " << parametersC.size() << " groups "	<< (*xp).size();
	}

	funObj.get_hessian(parametersC, resultC);
	if(configC.verbose) {
		Rcpp::Rcout  << "\n";
	}

	resultC.makeCompressed();

/*	Eigen::SparseMatrix<double> Hfull(
		resultC.selfadjointView<Eigen::Lower>()
		);*/
	Rcpp::S4 out = eigen_to_dgC(resultC);

	return(out);
}


#endif
