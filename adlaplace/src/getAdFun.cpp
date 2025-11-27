#include "adlaplace/adlaplace.hpp"
#include "adlaplace/adfun.hpp"
#include "adlaplace/lgamma.hpp"



template <class T>
struct nbSizeValues{
	T logNbSize;
	T nbSize;
	T lgammaNbSize;
	T sizeLogSize;

	explicit nbSizeValues<T>(const T logTheta) {
		logNbSize = -0.5*logTheta;
		nbSize = exp_any(logNbSize);
		lgammaNbSize = lgamma_any(nbSize); 
		sizeLogSize = nbSize * logNbSize;
	}

};

CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const Data& data,
	const Config& config,
	const size_t start, const size_t end
	){
	nbSizeValues nbSize(config.logTheta[config.logTheta.size()-1L]);
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


		const double expScale = eta > nbSize.logNbSize?-1.0:1.0;
		const double maxValue = eta > nbSize.logNbSize?CppAD::Value(eta):nbSize.logNbSize;
		CppAD::AD<double> diffEtaLogNbSize = eta - nbSize.logNbSize;
		CppAD::AD<double> logRplusMu = maxValue + CppAD::log1p(CppAD::exp(expScale*diffEtaLogNbSize));

		logDens += std::lgamma(data.y[Dobs] + nbSize.nbSize) 
		- std::lgamma(data.y[Dobs]  + 1.0)
		- nbSize.lgammaNbSize + nbSize.sizeLogSize 
		+ data.y[Dobs]  * eta
		- (nbSize.nbSize + data.y[Dobs]) * logRplusMu;

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
	const Data& data,
	const Config& config) {

	auto parameters = config.start_gamma;

	if(config.verbose) {
		Rcpp::Rcout << " adfun Nparams " << parameters.size() << " ";
	}

	std::vector<GroupPack> result(data.Ny+1L);

# pragma omp parallel
	{
		const size_t Nparams = parameters.size();
		CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
		for(size_t D=0;D<Nparams;D++) {
			ad_params[D] = parameters[D];
		}


# pragma omp for nowait
		for(size_t D=0;D<data.Ny;D++) {

			CppAD::Independent(ad_params);   
			auto resultHere = logDensObs(ad_params, data, config, D, D+1L);
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);
		}

# pragma omp single
		{
			CppAD::Independent(ad_params);
			auto resultHere = logDensRandom(ad_params, data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[data.Ny].fun = std::move(fun);
		}
	} // parallel

	if(config.group_sparsity.size()) {

		addSparsityToAdFun(result, config.group_sparsity);

	}

	return(result);
}




