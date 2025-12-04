
#include "adlaplace/adlaplace.hpp"
#include "adlaplace/lgamma.hpp"

#define COMPUTE_CONSTANTS

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

	double logNbSizeValue = to_double(logNbSize);


	const size_t startP = data.groups.p[Dgroup];
	const size_t endP = data.groups.p[Dgroup+1];

	for(size_t DI=startP;DI < endP; DI++) {
		const size_t Dobs = data.groups.i[DI];

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

// dummy funciton if there is no logDensExtra
#ifdef NOLOGDENSEXTRA
template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	) {
	CppAD::vector<Type> result(1);
	result[0] = 0.0;
	return(result);
}
#endif


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

	const size_t N=data.groups.i.size();
	for(size_t D=0; D <N;D++) {
		const size_t indexHere = data.groups.i[D];
		logDens1 += lgamma_any<Type>(data.y[indexHere] + nbSize); 
	}

	Type logDens2 = N*(sizeLogSize - lgammaNbSize);

#ifdef COMPUTE_CONSTANTS
// always constant	- std::lgamma(data.y[Dobs]  + 1.0)
	for(size_t D=0; D <N;D++) {
		const size_t indexHere = data.groups.i[D];
		logDens2 -= std::lgamma(data.y[indexHere] + 1.0); 
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


template <class Type>
CppAD::vector<CppAD::AD<double>> logDensRandom(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	){

	CppAD::vector<Type> logTheta(theta.size()), expTheta(theta.size());
	if(config.transform_theta) {
		for(size_t D=0;D<theta.size();++D) {
			logTheta[D] = theta[D];
			expTheta[D] = exp_any(theta[D]);
		}
	} else {
		for(size_t D=0;D<theta.size();++D) {
			logTheta[D] = log_any(theta[D]);
			expTheta[D] = theta[D];
		}
	}


	CppAD::vector<CppAD::AD<double>> result(1, 0.0);
	CppAD::AD<double> qpart = 0.0, qDet=0.0;

	CppAD::vector<CppAD::AD<double>> gammaScaled(data.Ngamma);


	if(config.verbose) {
		Rcpp::Rcout << "q, ngamma  " << data.Ngamma << " nmap " << data.map.size() << 
		" ntheta " << config.theta.size() << 
		" exp theta map 0 " <<  expTheta[data.map[0]] << " gamma0 " << gamma[0] << "\n";
	}


	for(size_t D=0;D<data.Ngamma;D++) {

		size_t mapHere = data.map[D];
		gammaScaled[D] = gamma[D] / expTheta[mapHere];
		qpart += 
			gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
		qDet += logTheta[mapHere];
	}
		qpart *= 0.5;

	if(config.verbose) {
		Rcpp::Rcout << "logDensRandom " << qpart << " det "<< qDet << "\n";
	}

	// Warning, no offdiag of Q implemented

	result[0] = qpart + qDet;
	return(result);
}



template CppAD::vector<CppAD::AD<double>>
logDensObs<double>(
    const CppAD::vector<CppAD::AD<double>>&,
    const CppAD::vector<double>&,
    const CppAD::vector<double>&,
    const Data&, const Config&, const size_t);

template CppAD::vector<CppAD::AD<double>>
logDensObs<CppAD::AD<double>>(
    const CppAD::vector<CppAD::AD<double>>&,
    const CppAD::vector<CppAD::AD<double>>&,
    const CppAD::vector<CppAD::AD<double>>&,
    const Data&, const Config&, const size_t);


template CppAD::vector<double>
logDensExtra<double>(
    const CppAD::vector<double>&,
    const Data&, const Config&);

template CppAD::vector<CppAD::AD<double>>
logDensExtra<CppAD::AD<double>>(
    const CppAD::vector<CppAD::AD<double>>&,
    const Data&, const Config&);


template CppAD::vector<CppAD::AD<double>>
logDensRandom<double>(
    const CppAD::vector<CppAD::AD<double>>&,
    const CppAD::vector<double>&,
    const Data&, const Config&);

template CppAD::vector<CppAD::AD<double>>
logDensRandom<CppAD::AD<double>>(
    const CppAD::vector<CppAD::AD<double>>&,
    const CppAD::vector<CppAD::AD<double>>&,
    const Data&, const Config&);



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
		Rcpp::Rcout << "getting ad function, inner " << inner << " ";
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