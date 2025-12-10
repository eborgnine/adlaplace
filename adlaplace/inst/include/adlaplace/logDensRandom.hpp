#ifndef LOGDENSRANDOM_HPP
#define LOGDENSRANDOM_HPP


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
		Rcpp::Rcout << "q adlaplace, ngamma  " << data.Ngamma << " nmap " << data.map.ncol() << 
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

	if(config.verbose) {
		Rcpp::Rcout << "logDensRandom " << qpart << " det "<< qDet << "\n";
	}

	// Warning, no offdiag of Q implemented

	result[0] = qpart + qDet;
	return(result);
}

// declare
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



#endif


