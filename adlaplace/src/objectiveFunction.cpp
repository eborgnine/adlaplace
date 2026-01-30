/*  Rcpp::depends(RcppEigen)]] */

#include "adlaplace/data.hpp"

// header for the lgamma function
#include "adlaplace/lgamma.hpp"

// use the standard logDensRandom from logDensRandom.hpp
#include "adlaplace/logDensRandom.hpp"

// returns log density of observations conditional on random effects

#define COMPUTE_CONSTANTS

CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<CppAD::AD<double>> &beta,
	const CppAD::vector<CppAD::AD<double>> &theta,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	){

	CppAD::AD<double> logDens1 = 0.0, logDens2 = 0.0;

	CppAD::AD<double> logTheta = config.transform_theta?
	theta[theta.size()-1]:CppAD::log(theta[theta.size()-1]);

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
			etaFixed += data.X.x[D] * beta[data.X.i[D]];
		}
		CppAD::AD<double> etaRandom = 0.0;
		for(size_t D =p0a;D < p1a; D++) {
			etaRandom += data.A.x[D] * gamma[data.A.i[D]];
		}

		const CppAD::AD<double> eta = etaRandom + etaFixed;

		// logNbSize is usually large
		const CppAD::AD<double> diffEtaLogNbSize = eta - logNbSize;
		const CppAD::AD<double> logRplusMu = logNbSize + CppAD::log1p(CppAD::exp(diffEtaLogNbSize));

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
	const CppAD::vector<CppAD::AD<double>> &theta,
	const Data& data,
	const Config& config
	) {

	CppAD::AD<double> logTheta = config.transform_theta?
	theta[theta.size()-1]:CppAD::log(theta[theta.size()-1]);

	CppAD::AD<double> logNbSize = -2*logTheta;
	CppAD::AD<double> nbSize = CppAD::exp(logNbSize);
	CppAD::AD<double>  lgammaNbSize = lgamma_ad(nbSize); 
	CppAD::AD<double>  sizeLogSize = nbSize * logNbSize;

	CppAD::AD<double>  logDens1=0.0;

	const size_t N=data.y.size();
	for(size_t D=0; D <N;D++) {
		logDens1 += lgamma_ad(data.y[D] + nbSize); 
	}

	CppAD::AD<double>  logDens2 = N*(sizeLogSize - lgammaNbSize);

// always constant	- std::lgamma(data.y[Dobs]  + 1.0)
#ifdef COMPUTE_CONSTANTS	
	double constants=0.0;
	for(size_t D=0; D <N;D++) {
		constants += std::lgamma(data.y[D] + 1.0); 
	}	
	logDens2 -= constants;
#endif
	CppAD::AD<double>  logDens = logDens1 + logDens2;
	if(config.verbose) {
		Rcpp::Rcout << "logDensExtra " << logDens << " " << logDens1 << " " << logDens2 << " nbSize " << nbSize << 
		" logTheta " << logTheta << " tr " << config.transform_theta << "\n";
	}
	CppAD::vector<CppAD::AD<double> > result(1);
	result[0] = logDens;
	return(result);
}

// ADfun and interfaces
#include"adlaplace/adlaplace_api.hpp"


