#include "adlaplace/adlaplace.hpp"
#include "adlaplace/adfun.hpp"
#include "adlaplace/foromp.hpp"
#include "adlaplace/innerOpt.hpp"

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

/*
These functions must be defined somewhere
std::vector<GroupPack> getAdFunInner(
  const Data& data,
  const Config& config);

std::vector<GroupPack> getAdFunOuter(
  const Data& data,
  const Config& config);
*/

//' @export
// [[Rcpp::export]]
SEXP getAdFun(
	Rcpp::List data, 
	Rcpp::List config,
	const bool inner=false)
{

	Data dataC(data);
	Config configC(config);
	if(configC.verbose) {
		Rcpp::Rcout << "getting ad function";
	}

	omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

	std::vector<GroupPack> adPack = inner?
		getAdFunInner(dataC, configC):getAdFunOuter(dataC, configC);
	auto* ptr = new std::vector<GroupPack>(std::move(adPack));
  Rcpp::XPtr<std::vector<GroupPack>> xp(ptr, /*deleteOnFinalizer=*/true);

	xp.attr("class") = "adpack_ptr";
	return xp;
}

//' @export
// [[Rcpp::export]]
double jointLogDens(
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

	AD_Func_Opt funObj(*xp, configC.hessianIPLower);

	const size_t Nparams = parameters.size();
	Eigen::VectorXd parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

	funObj.get_f(parametersC, result);

	return(result);

}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{
	Config configC(config);

	const size_t Nparams = parameters.size();
	Rcpp::NumericVector result(Nparams);
	if (adPack == R_NilValue) {
		return(Rcpp::NumericVector());
	}

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	AD_Func_Opt funObj(*xp, configC.hessianIPLower);


	Eigen::VectorXd parametersC(Nparams);
	Eigen::VectorXd resultC(Nparams);

	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	omp_set_num_threads(configC.num_threads);
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

//' @export
// [[Rcpp::export]]
Rcpp::S4 hessian(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	if (adPack == R_NilValue) {
		return(R_NilValue);
	}

	Config configC(config);

	const size_t Nparams = parameters.size();

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);

	AD_Func_Opt funObj(*xp, configC.hessianIPLower);



	Eigen::VectorXd parametersC(Nparams);
	
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	Eigen::SparseMatrix<double> resultC = funObj.Htemplate.cast<double>();

	omp_set_num_threads(configC.num_threads);
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

	Rcpp::S4 out = eigen_to_dgC(resultC);
	return(out);
}

//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
	Rcpp::NumericVector parameters, 
	Rcpp::List data,
	Rcpp::List control, 
	Rcpp::List config)
{

	Config configC(config);
	Data dataC(data);
	TrustControl ctrl(control); 

	using Tvec  = Eigen::VectorXd;

	const size_t Nparams = parameters.size();
	Tvec parametersC(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

	std::vector<GroupPack> adPack = getAdFunInner(dataC, configC);

	Rcpp::List res = inner_opt(parametersC, adPack, ctrl, configC);

	return(res);
}


//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt_adpack(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	const Rcpp::List control, 
	const Rcpp::List config)
{

	const Config configC(config);
	const TrustControl ctrl(control); 

	using Tvec  = Eigen::VectorXd;

	const size_t Nparams = parameters.size();
	Tvec parametersC(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);

	omp_set_num_threads(configC.num_threads);

	Rcpp::List res = inner_opt(parametersC, *xp, ctrl, configC);

	return(res);
}



