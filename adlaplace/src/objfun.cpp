#include "adlaplace/adlaplace.hpp"
#include "adlaplace/functions.hpp"
#include "adlaplace/foromp.hpp"

#include "configObj.hpp"


// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// from trustOptim
#include <common_R.hpp>
#include <CG-base.h>
#include <CG-quasi.h> 


CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const Data& data,
	const Config& config,
	const size_t start, const size_t end
	){

	CppAD::AD<double> logDens = 0.0;
	CppAD::vector<CppAD::AD<double>> result(1, 0.0);

	for(size_t Dobs=start;Dobs < end; Dobs++) {

		CppAD::AD<double>  eta = 0.0;
		const size_t p0x = data.X.p[Dobs];
		const size_t p1x = data.X.p[Dobs + 1];
		const size_t p0a = data.A.p[Dobs];
		const size_t p1a = data.A.p[Dobs + 1];
		for(size_t D =p0x;D < p1x; D++) {
			eta += data.X.x[D] * config.beta[data.X.i[D]];
		}
		for(size_t D =p0a;D < p1a; D++) {
			eta += data.A.x[D] * gamma[data.A.i[D]];
		}
		// poissondens = lambda^y exp(-lambda)/y!, lambda = exp(eta)
//		logDens += data.y[Dobs] * eta - CppAD::exp(eta);


		const double expScale = eta > config.logNbSize?-1.0:1.0;
		const double maxValue = eta > config.logNbSize?CppAD::Value(eta):config.logNbSize;
		CppAD::AD<double> diffEtaLogNbSize = eta - config.logNbSize;
		CppAD::AD<double> logRplusMu = maxValue + CppAD::log1p(CppAD::exp(expScale*diffEtaLogNbSize));

		logDens += std::lgamma(data.y[Dobs] + config.nbSize) 
		- std::lgamma(data.y[Dobs]  + 1.0)
		- config.lgammaNbSize + config.sizeLogSize 
		+ data.y[Dobs]  * eta
		- (config.nbSize + data.y[Dobs]) * logRplusMu;

	}
	result[0] = -logDens;
	return(result);
}

CppAD::vector<CppAD::AD<double>> logDensRandom(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const Data& data,
	const Config& config
	){

	CppAD::vector<CppAD::AD<double>> result(1, 0.0);
	CppAD::AD<double> qpart = 0.0;

	CppAD::vector<CppAD::AD<double>> gammaScaled(data.Ngamma);
	for(size_t D=0;D<data.Ngamma;D++) {
		size_t mapHere = data.map[D];
		gammaScaled[D] = gamma[D] / config.theta[mapHere];
		qpart += config.logTheta[mapHere] + 
		gammaScaled[D]*gammaScaled[D]*(0.5*data.Qdiag[D]);
	}

	// Warning, no offdiag of Q implemented

	result[0] = qpart;
	return(result);
}

std::vector<GroupPack> getAdFun(
	const CppAD::vector<double> & parameters,  
	const Data& data,
	const Config& config) {


	std::vector<GroupPack> result(data.Ny+1L);
	const size_t Nparams = parameters.size();

	CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		ad_params[D] = parameters[D];
	}

	for(size_t D=0;D<data.Ny;D++) {

		CppAD::Independent(ad_params);   
		auto resultHere = logDensObs(ad_params, data, config, D, D+1L);
		CppAD::ADFun<double> fun(ad_params, resultHere);
		result[D].fun = std::move(fun);
	}

	CppAD::Independent(ad_params);
	auto resultHere = logDensRandom(ad_params, data, config);   
	CppAD::ADFun<double> fun(ad_params, resultHere);
	result[data.Ny].fun = std::move(fun);

	return(result);
}


//' @export
// [[Rcpp::export]]
SEXP getAdFun(
	Rcpp::List data, 
	Rcpp::List config)
{

	Data dataC(data);
	Config configC(config);
	Rcpp::NumericVector parameters = config["start_gamma"];

	const size_t Nparams = parameters.size();

	CppAD::vector<double> parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	std::vector<GroupPack> adPack = getAdFun(parametersC, dataC, configC);

	if(config.containsElementNamed("group")) {
		if(configC.verbose) {
			Rcpp::Rcout << "adding sparsity patterns\n";	 
		}
		Rcpp::List group_sparsity_list = config["group"];
		addSparsityToAdFun(adPack, group_sparsity_list);
	}

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
	std::vector<GroupPack>& adPackC = *xp;

	AD_Func_Opt funObj(adPackC);


	const size_t Nparams = parameters.size();
	Eigen::VectorXd parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}


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
		return(result);
	}

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
//	std::vector<GroupPack>& adPackC = *xp;

	AD_Func_Opt funObj(*xp);


	Eigen::VectorXd parametersC(Nparams);
	Eigen::VectorXd resultC(Nparams);

	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

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
// or equivalently:
// Eigen::SparseMatrix<double> resultC(funObj.Htemplate.cast<double>());

	funObj.get_hess(parametersC, resultC);

resultC.makeCompressed();

Rcpp::S4 out = eigen_to_dgC(resultC);
	return(out);

}


//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List control, 
	Rcpp::List config)
{

		Config configC(config);

	using Tvec   = Eigen::VectorXd;
	using THess   = Eigen::MatrixXd;
	using TPreLLt = Eigen::LLT<THess>;

	const size_t Nparams = parameters.size();
	Tvec parametersC(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}



	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	AD_Func_Opt funObj(*xp, configC.hessianIPLower);

	TrustControl ctrl(control); 

	Trust_CG_Optimizer<Tvec, AD_Func_Opt, THess, TPreLLt> opt(
		funObj, 
		parametersC,
		ctrl.rad,
		ctrl.min_rad,
		ctrl.tol,
		ctrl.prec,
		ctrl.report_freq,
		ctrl.report_level,
		ctrl.header_freq,
		ctrl.report_precision,
		ctrl.maxit,
		ctrl.contract_factor,
		ctrl.expand_factor,
		ctrl.contract_threshold,
		ctrl.expand_threshold_rad,
		ctrl.expand_threshold_ap,
		ctrl.function_scale_factor,
		ctrl.precond_refresh_freq,
		ctrl.precond_ID,
		ctrl.quasi_newton_method,
		ctrl.trust_iter
		);


  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

	opt.run();

    // collect results and return

	Tvec P(Nparams);
	Tvec grad(Nparams);

	double fval, radius;
	int iterations;
	MB_Status status;

	status = opt.get_current_state(P, fval, grad, iterations, radius);

	Rcpp::List res = Rcpp::List::create(
		Rcpp::Named("fval") = Rcpp::wrap(fval),
		Rcpp::Named("solution") = Rcpp::wrap(P),
		Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		Rcpp::Named("method") = Rcpp::wrap("Sparse")
		);

	return(res);


}
