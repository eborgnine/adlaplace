#include "adlaplace/adlaplace.hpp"

//#define DEBUG

// some constants
// 0.5 * log(2*pi)
const double ONEHALFLOGTWOPI = 0.9189385332046727417803297364056176398613974736377834128171515404;
//  sqrt(2), std::numbers::sqrt2 isn't always available
const double SQRTTWO = 1.414213562373095048801688724209698078569671875376948073176679737990732478462107;
static const std::string JAC_COLOR = "cppad";  
static const std::string HESS_COLOR = "cppad.symmetric";


// use logDensRandom from logDensRandom.hpp
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

	CppAD::AD<double> result = 0.0;
	CppAD::vector<CppAD::AD<double>> resultV(1, 0.0);

	// last two elements omega and alpha
	Type omega = config.transform_theta ? 
	exp_any(theta[theta.size()-2L]) : 
	theta[theta.size()-2L];

	Type omegaTimesSqrtTwo = omega * Type(SQRTTWO);

	Type alpha = theta[theta.size()-1L];

	const bool have_groups = config.groups.ncol() > 0;
	const size_t startP = have_groups?config.groups.p[Dgroup]:Dgroup;
	const size_t endP = have_groups?config.groups.p[Dgroup+1]:Dgroup+1;

	for(size_t DI=startP;DI < endP; DI++) {
		const size_t Dobs = have_groups?config.groups.i[DI]:DI;

		CppAD::AD<double>  etaRandom = 0.0;
		Type  etaFixed = 0.0;
		const size_t p1x = data.X.p[Dobs + 1];
		for(size_t D =data.X.p[Dobs];D < p1x; D++) {
			etaFixed += data.X.x[D] * beta[data.X.i[D]];
		}
		const size_t p1a = data.A.p[Dobs + 1];
		for(size_t D =data.A.p[Dobs];D < p1a; D++) {
			etaRandom += data.A.x[D] * gamma[data.A.i[D]];
		}

		const CppAD::AD<double> eta = etaFixed + etaRandom;
		const CppAD::AD<double> z = (data.y[Dobs] - eta) / omegaTimesSqrtTwo;
		const CppAD::AD<double> erfcArg = -alpha * z;
		const CppAD::AD<double> erfcVal = CppAD::erfc(erfcArg);
		const CppAD::AD<double> toAdd = z*z - CppAD::log(erfcVal);

#ifdef DEBUG
		Rcpp::Rcout << "Dobs " << Dobs << " eta " << eta << " z " << z << " erfcArg " << erfcArg <<
		" erfcVal " << erfcVal << " toAdd " << toAdd << "\n";
#endif
		result += toAdd; 
	}
#ifdef DEBUG
	Rcpp::Rcout << "Dgroup " << Dgroup << "result " << result << "\n";
#endif

	resultV[0] = result;
	return(resultV);
}



template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	) {

	Type logOmega = config.transform_theta?
	theta[theta.size()-2L]:
	log_any(theta[theta.size()-2L]);

	const size_t N=data.y.size();

	Type logDens = Type(N)*logOmega + Type(N*ONEHALFLOGTWOPI);

	if(config.verbose) {
		Rcpp::Rcout << "logDensExtra " << logDens << " " <<
		" logOmega " << logOmega << " tr " << config.transform_theta << "\n";
	}
	CppAD::vector<Type> result(1);
	result[0] = logDens;
	return(result);
}


// declare
#include"adlaplace/declaration_macros.hpp"


