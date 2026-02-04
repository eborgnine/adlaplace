#include "adlaplace/utils.hpp"

// header for the lgamma function
#include "adlaplace/lgamma.hpp"


// returns log density of observations conditional on random effects

//#define DEBUG can't multi-thread with debug enabled

// use the standard logDensRandom from logDensRandom.hpp
#include "adlaplace/logDensRandom.hpp"


CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& x,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	){

	CppAD::AD<double> logDens1 = 0.0, logDens2 = 0.0;

	CppAD::AD<double> thetaIn = x[config.theta_end-1];
	CppAD::AD<double> logTheta = config.transform_theta?thetaIn:CppAD::log(thetaIn);

	CppAD::AD<double> logNbSize = -2*logTheta;
	CppAD::AD<double> nbSize = CppAD::exp(logNbSize);

	const bool have_groups = config.groups.ncol() > 0;
	const size_t startP = have_groups?config.groups.p[Dgroup]:Dgroup;
	const size_t endP = have_groups?config.groups.p[Dgroup+1]:Dgroup+1;

	for(size_t DI=startP;DI < endP; DI++) {
		const size_t Dobs = have_groups?config.groups.i[DI]:DI;

		const size_t p0x = data.X.p[Dobs];
		const size_t p1x = data.X.p[Dobs + 1];
		const size_t p0a = data.A.p[Dobs];
		const size_t p1a = data.A.p[Dobs + 1];
		CppAD::AD<double>  etaFixed = 0.0;
		for(size_t D =p0x;D < p1x; D++) {
			etaFixed += data.X.x[D] * x[config.beta_begin + data.X.i[D]];
		}
		CppAD::AD<double> etaRandom = 0.0;
		for(size_t D =p0a;D < p1a; D++) {
			etaRandom += data.A.x[D] * x[config.gamma_begin + data.A.i[D]];
		}

		const CppAD::AD<double> eta = etaRandom + etaFixed;

		const CppAD::AD<double> diff = eta - logNbSize;
//		const CppAD::AD<double> logRplusMu = logNbSize + CppAD::log1p(CppAD::exp(diffEtaLogNbSize));
// softplus(diff) = log(1 + exp(diff)) computed stably
		const CppAD::AD<double> softplus = CppAD::CondExpGt(
			diff, CppAD::AD<double>(0.0),
			diff + CppAD::log1p(CppAD::exp(-diff)),   
			CppAD::log1p(CppAD::exp(diff))     
			);

		const CppAD::AD<double> logRplusMu = logNbSize + softplus;


		logDens1 += data.y[Dobs] * eta;
		logDens2 += - logRplusMu * (data.y[Dobs] + nbSize);
// this bit is done in logDensExtra
/*		logDens +=
			+ lgamma(data.y[Dobs] + nbSize.nbSize) 
			- std::lgamma(data.y[Dobs]  + 1.0)
			- nbSize.lgammaNbSize + nbSize.sizeLogSize 
*/
	}
	CppAD::vector<CppAD::AD<double>> result(1);
	result[0] = logDens1 + logDens2;
	return(result);
}


CppAD::vector<CppAD::AD<double>> logDensExtra(
	const CppAD::vector<CppAD::AD<double>> &x,
	const Data& data,
	const Config& config
	) {

	CppAD::AD<double> thetaIn = x[config.theta_end-1];
	CppAD::AD<double> logTheta = config.transform_theta?thetaIn:CppAD::log(thetaIn);

	CppAD::AD<double> logNbSize = -2*logTheta;
	CppAD::AD<double> nbSize = CppAD::exp(logNbSize);
	CppAD::AD<double>  lgammaNbSize = lgamma_ad(nbSize); 
	CppAD::AD<double>  sizeLogSize = nbSize * logNbSize;

	CppAD::AD<double>  logDens1=0.0;

	for(size_t D=0; D < data.Ny; D++) {
		logDens1 += lgamma_ad(data.y[D] + nbSize); 
	}

	CppAD::AD<double>  logDens2 = data.Ny*(sizeLogSize - lgammaNbSize);

// always constant	- std::lgamma(data.y[Dobs]  + 1.0)
	double constants=0.0;
	for(size_t D=0; D <data.Ny; D++) {
		constants += std::lgamma(data.y[D] + 1.0); 
	}	
	logDens2 -= constants;

	CppAD::AD<double>  logDens = logDens1 + logDens2;
#ifdef DEBUG
	if(config.verbose) {
		Rcpp::Rcout << "logDensExtra " << logDens << " " << logDens1 << " " << logDens2 << " nbSize " << nbSize << 
		" logTheta " << logTheta << " tr " << config.transform_theta << "\n";
	}
#endif	
	CppAD::vector<CppAD::AD<double> > result(1);
	result[0] = logDens;
	return(result);
}

// ADfun and interfaces
#include"adlaplace/adlaplace_create.hpp"



