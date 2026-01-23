/*  Rcpp::depends(RcppEigen)]] */

#include "adlaplace/adlaplace.hpp"

// header for the lgamma function
#include "adlaplace/lgamma.hpp"

// use the standard logDensRandom from logDensRandom.hpp
#include "adlaplace/logDensRandom.hpp"

// returns minus log density
template <class Type>
CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &beta,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	){

	CppAD::AD<double> logDens1 = 0.0, logDens2 = 0.0;
	CppAD::vector<CppAD::AD<double>> result(1, 0.0);

	Type logTheta = config.transform_theta?
		theta[theta.size()-1]:log_any<Type>(theta[theta.size()-1]);

	Type logNbSize = -2*logTheta;
	Type nbSize = exp_any(logNbSize);

//	double logNbSizeValue = to_double(logNbSize);

	const bool have_groups = config.groups.ncol() > 0;
	const size_t startP = have_groups?config.groups.p[Dgroup]:Dgroup;
	const size_t endP = have_groups?config.groups.p[Dgroup+1]:Dgroup+1;

	for(size_t DI=startP;DI < endP; DI++) {
		const size_t Dobs = have_groups?config.groups.i[DI]:DI;

		CppAD::AD<double>  etaRandom = 0.0;
		Type  etaFixed = 0.0;
		const size_t p0x = data.X.p[Dobs];
		const size_t p1x = data.X.p[Dobs + 1];
		const size_t p0a = data.A.p[Dobs];
		const size_t p1a = data.A.p[Dobs + 1];
		for(size_t D =p0x;D < p1x; D++) {
			etaFixed += data.X.x[D] * beta[data.X.i[D]];
		}
		for(size_t D =p0a;D < p1a; D++) {
			etaRandom += data.A.x[D] * gamma[data.A.i[D]];
		}

		const CppAD::AD<double> eta = etaFixed + etaRandom;

//		const bool etaIsBigger = eta > logNbSize;
//		const double expScale = etaIsBigger?-1.0:1.0;
//		const double max_value = etaIsBigger?CppAD::Value(eta):logNbSizeValue;
//		const CppAD::AD<double> logRplusMu = max_value + CppAD::log1p(CppAD::exp(expScale*diffEtaLogNbSize));
		// logNbSize is usually large
		const CppAD::AD<double> diffEtaLogNbSize = eta - logNbSize;
		const CppAD::AD<double> logRplusMu = logNbSize + CppAD::log1p(CppAD::exp(diffEtaLogNbSize));

		logDens1 += data.y[Dobs]  * eta;
		logDens2 += - logRplusMu * (data.y[Dobs] + nbSize);
// this bit is done in logDensExtra
/*		logDens +=
			+ lgamma(data.y[Dobs] + nbSize.nbSize) 
			- std::lgamma(data.y[Dobs]  + 1.0)
			- nbSize.lgammaNbSize + nbSize.sizeLogSize 
*/

		}
//Rcpp::Rcout << "ld1 " << logDens1 << " ld2 " << logDens2 << "\n";		

	result[0] = -logDens1 - logDens2;
	return(result);
}



template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	) {

	Type logTheta = config.transform_theta?
		theta[theta.size()-1]:log_any<Type>(theta[theta.size()-1]);

	Type logNbSize = -2*logTheta;
	Type nbSize = exp_any<Type>(logNbSize);
	Type lgammaNbSize = lgamma_any<Type>(nbSize); 
	Type sizeLogSize = nbSize * logNbSize;


	Type logDens1=0.0;

	const size_t N=data.y.size();
	for(size_t D=0; D <N;D++) {
		logDens1 += lgamma_any<Type>(data.y[D] + nbSize); 
	}

	Type logDens2 = N*(sizeLogSize - lgammaNbSize);

#ifdef COMPUTE_CONSTANTS
// always constant	- std::lgamma(data.y[Dobs]  + 1.0)
	for(size_t D=0; D <N;D++) {
		logDens2 -= std::lgamma(data.y[D] + 1.0); 
	}	
#endif
	Type logDens = logDens1 + logDens2;
	if(config.verbose) {
		Rcpp::Rcout << "logDensExtra " << logDens << " " << logDens1 << " " << logDens2 << " nbSize " << nbSize << 
			" logTheta " << logTheta << " tr " << config.transform_theta << "\n";
	}
	CppAD::vector<Type> result(1);
	result[0] = -logDens;
	return(result);
}

// declare
#include"adlaplace/declaration_macros.hpp"


