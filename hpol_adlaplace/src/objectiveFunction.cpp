#define DEBUG

#include "adlaplace/adlaplace_base.hpp"
#include "adlaplace/lgamma.hpp"


//#define COMPUTE_CONSTANTS



CppAD::AD<double> stable_logsumexp(const CppAD::vector<CppAD::AD<double>>& eta) {
    // compute log(sum(exp(eta)))

    // find index of maximum
    size_t max_idx = 0;
    for (size_t Deta = 1; Deta < eta.size(); ++Deta) {
        if (eta[Deta] > eta[max_idx]) {
            max_idx = Deta;
        }
    }

    // get max value, sever tape if AD
    double max_value = CppAD::Value(eta[max_idx]);

    // sum exp differences
    CppAD::AD<double> sumexp = 0.0;
    CppAD::AD<double> max_valueT = max_value;

    for (size_t Deta = 0; Deta < eta.size(); ++Deta) {
        sumexp += CppAD::exp(eta[Deta] - max_valueT);
    }

    CppAD::AD<double> logSum = CppAD::log(sumexp);
    CppAD::AD<double> result = logSum + max_valueT;

    return(result);
}


// Compute eta for one stratum,  
template <class BetaT>
CppAD::vector<CppAD::AD<double>> compute_eta_for_stratum(size_t Dstrata,
                             const Data& data,
                             const CppAD::vector<CppAD::AD<double>>& gamma,
                             const CppAD::vector<BetaT>& beta)
{

  const size_t startHere = data.elgm_matrix.p[Dstrata], endHere = data.elgm_matrix.p[Dstrata+1];
  const size_t NinStrata = endHere - startHere;

  CppAD::vector<CppAD::AD<double>> etaHere(NinStrata);


  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const int Deta = data.elgm_matrix.i[k];

    CppAD::AD<double> accGamma = 0;
    BetaT accBeta = 0;

    // X contribution
    const int p0x = data.X.p[Deta];
    const int p1x = data.X.p[Deta + 1];
    for (int t = p0x; t < p1x; ++t) {
      accBeta += data.X.x[t] * beta[data.X.i[t]];
    }

    // A contribution
    const int p0a = data.A.p[Deta];
    const int p1a = data.A.p[Deta + 1];
    for (int t = p0a; t < p1a; ++t) {
      accGamma += data.A.x[t] * gamma[data.A.i[t]];
    }

    etaHere[j] = accBeta + accGamma;
  }
  return(etaHere);
}

// Accumulate contribution over one stratum.
template <class ThetaT>
CppAD::AD<double> accumulate_contrib_for_stratum(size_t Dstrata,
                               const Data& data,
                               const CppAD::vector<CppAD::AD<double>>& etaHere, // length = NinStrata
							   const ThetaT &logNuSq)
{
  const size_t startHere = data.elgm_matrix.p[Dstrata];
  const size_t endHere   = data.elgm_matrix.p[Dstrata + 1];
  const size_t NinStrata = endHere - startHere;


  CppAD::AD<double> contrib = 0;

  if (NinStrata == 0) {
    return(contrib);
  }

  // normalize to get probabilities
  CppAD::AD<double> etaLogSum = stable_logsumexp(etaHere);

  for (size_t j = 0, k = startHere; j < NinStrata; ++j, ++k) {
    const size_t idx = static_cast<size_t>(data.elgm_matrix.i[k]);
    const double yhere = data.y[idx];

//    CppAD::AD<double> yOut = yhere;

    const CppAD::AD<double> etaMinusLogSumMu = etaHere[j] - etaLogSum;
    const CppAD::AD<double> muBarDivNuSq   = CppAD::exp(etaMinusLogSumMu - logNuSq);

      contrib += lgamma_ad(muBarDivNuSq + yhere) - lgamma_ad(muBarDivNuSq);    
  }

  return contrib;
}


template <class Type>
CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &beta,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	){

	const bool have_groups = config.groups.ncol() > 0;
	const size_t startP = have_groups?config.groups.p[Dgroup]:Dgroup;
	const size_t endP = have_groups?config.groups.p[Dgroup+1]:Dgroup+1;

	const Type lastTheta = theta[theta.size()-1];
	const Type logNuSq = config.transform_theta?2*lastTheta:2*exp_any(lastTheta);

	CppAD::AD<double> result1=0.0;
	CppAD::vector<CppAD::AD<double>> result(1);

	for(size_t DstrataI=startP;DstrataI < endP; DstrataI++) {

		const size_t Dstrata = have_groups?config.groups.i[DstrataI]:DstrataI;

  		auto etaHere = compute_eta_for_stratum(
    		Dstrata, data, gamma, beta);

  		auto contrib = accumulate_contrib_for_stratum(
    		Dstrata, data, etaHere, logNuSq);

//  		Rcpp::Rcout << "D " << Dstrata << " e " << etaHere[0] << " c " << contrib << " ";
		result1 += contrib;
	}
	result[0] = -result1;
	return(result);
}


template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	) {

	const Type lastTheta = theta[theta.size()-1];
	const Type logNuSq = config.transform_theta?2*lastTheta:2*exp_any(lastTheta);
    const Type oneOverNuSq  = CppAD::exp(-logNuSq);

    
    const size_t Nstrata = data.elgm_matrix.ncol();

    Type contribLgamma1overNuSqPlusSumYil=0;
#ifdef COMPUTE_CONSTANTS    
    Type contribLgamma1pSumYil=0;
    Type contribLgammaYp1=0;
#endif


	if(config.verbose) {
		Rcpp::Rcout << "extra Nstrata " << Nstrata << " lastTheta " <<
		lastTheta << " oneOverNuSq " << oneOverNuSq;
	}

    for(size_t Dstrata = 0;Dstrata<Nstrata;Dstrata++) {
    	int sumYhere = 0;
    	const size_t endHere = data.elgm_matrix.p[Dstrata+1];
    	for(size_t DobsI = data.elgm_matrix.p[Dstrata];DobsI < endHere;DobsI++) {
    		const size_t Dobs = data.elgm_matrix.i[DobsI];
    		sumYhere += data.y[Dobs];
#ifdef COMPUTE_CONSTANTS    		
    		contribLgammaYp1 += std::lgamma(1 + data.y[Dobs]);
#endif    		
    	}
    	contribLgamma1overNuSqPlusSumYil += lgamma_any(oneOverNuSq + sumYhere);
#ifdef COMPUTE_CONSTANTS  
    	contribLgamma1pSumYil += std::lgamma(1 + sumYhere);
#endif    	
    }


	if(config.verbose) {
		Rcpp::Rcout << " contribLgamma1overNuSqPlusSumYil " << contribLgamma1overNuSqPlusSumYil <<
		" lgamma_any(oneOverNuSq) " << lgamma_any(oneOverNuSq);
	}


   	Type contrib = Nstrata * lgamma_any(oneOverNuSq) + contribLgamma1overNuSqPlusSumYil;
 
#ifdef COMPUTE_CONSTANTS  
    contrib += contribLgammaYp1 + contribLgamma1pSumYil;
#endif    	

	if(config.verbose) {
		Rcpp::Rcout << " contrib " << contrib << "\n";
	}

	CppAD::vector<Type> result(1);
	result[0] = contrib;
	return(result);
}
 
template <class Type>
CppAD::vector<CppAD::AD<double>> logDensRandom(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	){

	const size_t Ngamma = gamma.size();
	const size_t Nmap = data.map.size();
	const size_t startGamma = Ngamma - Nmap;

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



	CppAD::vector<CppAD::AD<double>> gammaScaled(Ngamma);


	if(config.verbose) {
		Rcpp::Rcout << "q, Ngamma  " << Ngamma << " nmap " << data.map.size() << 
		" ntheta " << config.theta.size() << 
		" exp theta map 0 " <<  expTheta[data.map[0]] << " gamma0 " << gamma[0] << "\n";
	}


	for(size_t D=0;D<startGamma;D++) {
		gammaScaled[D] = gamma[D];
		qpart += 
			gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
	}
	for(size_t D=startGamma,Dmap=0;D<Ngamma;D++,Dmap++) {
		Rcpp::Rcout << D;
		size_t mapHere = data.map[Dmap];
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


// declare

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


