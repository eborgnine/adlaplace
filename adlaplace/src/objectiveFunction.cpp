
#include "adlaplace/adlaplace.hpp"
#include "adlaplace/lgamma.hpp"


template <class Type>
CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &beta,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	){

	CppAD::AD<double> logDens = 0.0;
	CppAD::vector<CppAD::AD<double>> result(1, 0.0);

	Type logNbSizeTheta = config.transform_theta?
		theta[theta.size()]:log_any(theta[theta.size()]);

	Type logNbSize = -0.5*logNbSizeTheta;
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

		const bool etaIsBigger = eta > logNbSize;
		const double expScale = etaIsBigger?-1.0:1.0;
		const double max_value = etaIsBigger?CppAD::Value(eta):logNbSizeValue;

		const CppAD::AD<double> diffEtaLogNbSize = eta - logNbSize;
		const CppAD::AD<double> logRplusMu = max_value + CppAD::log1p(CppAD::exp(expScale*diffEtaLogNbSize));

		logDens += data.y[Dobs]  * eta
			 + (data.y[Dobs] - nbSize) * logRplusMu;
// this bit is done in logDensExtra
/*		logDens +=
			+ lgamma(data.y[Dobs] + nbSize.nbSize) 
			- std::lgamma(data.y[Dobs]  + 1.0)
			- nbSize.lgammaNbSize + nbSize.sizeLogSize 
*/

		}




	result[0] = -logDens;
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
		theta[theta.size()]:log_any(theta[theta.size()]);

	Type logNbSize = -0.5*logTheta;
	Type nbSize = exp_any(logNbSize);
	Type lgammaNbSize = lgamma_any(nbSize); 
	Type sizeLogSize = nbSize * logNbSize;


	Type logDens=0.0;

	const size_t N=data.groups.i.size();
	for(size_t D=0; D <N;D++) {
		const size_t indexHere = data.groups.i[D];
		logDens += lgamma_any(data.y[indexHere] + nbSize); 
	}
	logDens += N*(- lgammaNbSize + sizeLogSize);

// always constant	- std::lgamma(data.y[Dobs]  + 1.0)

	if(config.verbose) {
		Rcpp::Rcout << "logDensExtra " << logDens << "\n";
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
	CppAD::AD<double> qpart = 0.0;

	CppAD::vector<CppAD::AD<double>> gammaScaled(data.Ngamma);


	if(config.verbose) {
		Rcpp::Rcout << "q, ngamma  " << data.Ngamma << " nmap " << data.map.size() << " ntheta " << config.theta.size() << "\n";
	}


	for(size_t D=0;D<data.Ngamma;D++) {


		size_t mapHere = data.map[D];
		gammaScaled[D] = gamma[D] / theta[mapHere];
		qpart += logTheta[mapHere] + 
			gammaScaled[D]*gammaScaled[D]*(0.5*data.Qdiag[D]);

	}

	if(config.verbose) {
		Rcpp::Rcout << "logDensRandom " << qpart << "\n";
	}

	// Warning, no offdiag of Q implemented

	result[0] = qpart;
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