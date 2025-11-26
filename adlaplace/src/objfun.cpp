#include "adlaplace/adlaplace.hpp"
#include "adlaplace/trustOptimUtils.hpp"
#include "adlaplace/lgamma.hpp"
#include "adlaplace/adfun.hpp"

#include "configObj.hpp"

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// from trustOptim
#include <common_R.hpp>
#include <CG-base.h>
#include <CG-quasi.h> 

#define DEBUG

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
	Rcpp::NumericVector parameters, 
	Rcpp::List data, 
	Rcpp::List config)
{

	Data dataC(data);
	Config configC(config);

	const size_t Nparams = parameters.size();

	CppAD::vector<double> parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	std::vector<GroupPack> adpack = getAdFun(parametersC, dataC, configC);
	auto* ptr = new std::vector<GroupPack>(std::move(adpack));
  Rcpp::XPtr<std::vector<GroupPack>> xp(ptr, /*deleteOnFinalizer=*/true);

	xp.attr("class") = "adpack_ptr";
	return xp;
}

void grad(
	const CppAD::vector<double> &x,
	double &f,
	std::vector<double> &g, 
	std::vector<GroupPack>& adpack
	) {

	const std::size_t n = static_cast<std::size_t>(x.size());

	double fOut = 0.0;
	std::vector<double> gOut(n,0.0);
	CppAD::vector<double> w(1);
	w[0] = 1.0;

	for(size_t D=0;D<adpack.size();D++) {
		CppAD::vector<double> y = adpack[D].fun.Forward(0, x);
		fOut += y[0];
		CppAD::vector<double> gradHere = adpack[D].fun.Reverse(1, w);
		for(size_t Dpar=0;Dpar<n;Dpar++) {
			gOut[Dpar] += gradHere[Dpar];
		}
	}

	f = fOut;
	for(size_t Dpar=0;Dpar<n;Dpar++) {
		g[Dpar] = gOut[Dpar];
	}
}

double jointLogDens(
	const CppAD::vector<double> &x,
	std::vector<GroupPack>& adpack
	) {

	double fOut = 0.0;

	for(size_t D=0;D<adpack.size();D++) {
		CppAD::vector<double> y = adpack[D].fun.Forward(0, x);
		fOut += y[0];
	}
	return(fOut);
}

struct AD_Func_Opt {
	std::vector<GroupPack> &tape;

	AD_Func_Opt(std::vector<GroupPack> &tape_) : tape(tape_) {}

    // f and g together
    template <class DerivedX, class DerivedG>
	void get_fdf(const Eigen::MatrixBase<DerivedX> &x,
		double &f,
		Eigen::MatrixBase<DerivedG> &g) {
		const std::size_t n = static_cast<std::size_t>(x.size());
		std::vector<double> gOut(n,0.0);
		CppAD::vector<double> xp(n);
		for (std::size_t i = 0; i < n; ++i) {
			xp[i] = x[i];
		}

		grad(xp, f, gOut, tape);
		for (std::size_t i = 0; i < n; ++i) {
			g[i] = gOut[i];
		}
	}

    // f only
    template <class DerivedX>
	void get_f(const Eigen::MatrixBase<DerivedX> &x,
		double &f) {
		const std::size_t n = static_cast<std::size_t>(x.size());
		CppAD::vector<double> xp(n);
		for (std::size_t i = 0; i < n; ++i) {
			xp[i] = x[i];
		}
		f = jointLogDens(xp, tape);

	}

    // g only
    template <class DerivedX, class DerivedG>
	void get_df(const Eigen::MatrixBase<DerivedX> &x,
		Eigen::MatrixBase<DerivedG> &g) {
		const std::size_t n = static_cast<std::size_t>(x.size());
		std::vector<double> gOut(n,0.0);
		CppAD::vector<double> xp(n);
		for (std::size_t i = 0; i < n; ++i) {
			xp[i] = x[i];
		}
		double f;
		grad(xp, f, gOut, tape);
		for (std::size_t i = 0; i < n; ++i) {
			g[i] = gOut[i];
		}
	}
};


//' @export
// [[Rcpp::export]]
double jointLogDens(
	Rcpp::NumericVector parameters, 
	SEXP adPack)
{
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
	SEXP adPack)
{
	const size_t Nparams = parameters.size();
	Rcpp::NumericVector result(Nparams);
	if (adPack == R_NilValue) {
		return(result);
	}

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	std::vector<GroupPack>& adPackC = *xp;

	AD_Func_Opt funObj(adPackC);


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
Rcpp::List inner_opt(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List control)
{

	using Tvec   = Eigen::VectorXd;
	using THess   = Eigen::MatrixXd;
	using TPreLLt = Eigen::LLT<THess>;

	const size_t Nparams = parameters.size();
	Tvec parametersC(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	std::vector<GroupPack>& adPackC = *xp;

	AD_Func_Opt funObj(adPackC);

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
