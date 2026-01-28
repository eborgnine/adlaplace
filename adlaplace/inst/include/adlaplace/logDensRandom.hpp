#ifndef LOGDENSRANDOM_HPP
#define LOGDENSRANDOM_HPP

const double ONEHALFLOGTWOPIRANDOM = 0.9189385332046727417803297364056176398613974736377834128171515404;

/*
	the standard (negative) log density for random effects
	include with #include "adlaplace/logDensRandom.hpp"
*/

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
		Rcpp::Rcout << "q ngamma  " << data.Ngamma << " map.ncol " << data.map.ncol() << 
		" map@i.size " << data.map.i.size() << 
		" ntheta " << config.theta.size() << 
		" exp theta map 0 " <<  expTheta[data.map.i[0]] << " gamma0 " << gamma[0] << ".\n";

	}


	for(size_t D=0;D<data.Ngamma;D++) {
		gammaScaled[D] = gamma[D];
	}

	for(size_t Dtheta=0;Dtheta<data.Nmap;Dtheta++) {
		const size_t mapStart = data.map.p[Dtheta];
		const size_t mapEnd = data.map.p[Dtheta+1];	
		const size_t Nhere = mapEnd - mapStart;

		if(Nhere) {
			qDet += logTheta[Dtheta] * Nhere;
		}

		for(size_t DgammaI=mapStart;DgammaI<mapEnd;DgammaI++) {
			const size_t Dgamma = data.map.i[DgammaI];
			gammaScaled[Dgamma] /= expTheta[Dtheta];
		}
	}

	for(size_t D=0;D<data.Ngamma;D++) {
		qpart += 
			gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
	}
	qpart *= 0.5;

	// Warning, no offdiag of Q implemented

	result[0] = qpart + qDet + CppAD::AD<double>(data.Ngamma * ONEHALFLOGTWOPIRANDOM);

	if(config.verbose) {
		Rcpp::Rcout << "logDensRandom " << result[0] << " qpart " << qpart << " det "<< qDet << "\n";
	}



	return(result);
}

#endif


