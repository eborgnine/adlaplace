#ifndef LOGDENSRANDOM_HPP
#define LOGDENSRANDOM_HPP

#include "adlaplace/constants.hpp"

/*
	the standard log density for random effects
	include with #include "adlaplace/logDensRandom.hpp"
*/

CppAD::vector<CppAD::AD<double>> logDensRandom(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<CppAD::AD<double>> &theta,
	const Data& data,
	const Config& config
	){

	CppAD::vector<CppAD::AD<double>> logTheta(theta.size()), expTheta(theta.size());
	if(config.transform_theta) {
		for(size_t D=0;D<theta.size();++D) {
			logTheta[D] = theta[D];
			expTheta[D] = CppAD::exp(theta[D]);
		}
	} else {
		for(size_t D=0;D<theta.size();++D) {
			logTheta[D] = CppAD::log(theta[D]);
			expTheta[D] = theta[D];
		}
	}

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
#ifdef COMPUTE_CONSTANTS	
	qDet += CppAD::AD<double>(data.Ngamma * ONEHALFLOGTWOPI)
#endif	

	CppAD::vector<CppAD::AD<double>> result(1, 0.0);
	result[0] = - qpart - qDet;

	if(config.verbose) {
		Rcpp::Rcout << "logDensRandom " << result[0] << " qpart " << qpart << " det "<< qDet << "\n";
	}
	return(result);
}

#endif


