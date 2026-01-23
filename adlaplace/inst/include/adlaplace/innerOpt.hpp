#ifndef INNER_OPT_HPP
#define INNER_OPT_HPP


inline Rcpp::List inner_opt(
	Eigen::VectorXd& parameters, 
	std::vector<GroupPack>& adPack,
	const TrustControl& ctrl, 
	const Config& config,
	std::vector<GroupPack>* adPackFull = nullptr)
{
	using Tvec   = Eigen::VectorXd;
	using THess   = Eigen::SparseMatrix<double>; 
	using TPreLLt = Eigen::SimplicialLLT<THess>;

	if(config.verbose) {
		Rcpp::Rcout << "getting ad fun object ";
	}


	AD_Func_Opt funObj(adPack, config.hessianIPLower_inner);

	if(config.verbose) {
		Rcpp::Rcout << "starting opt ";
	}

	Trust_CG_Sparse<Tvec, AD_Func_Opt, THess, TPreLLt> opt(
		funObj, 
		parameters,
		ctrl.rad,
		ctrl.min_rad,
		ctrl.tol,
		ctrl.prec,
		ctrl.report_freq,
		ctrl.report_level,
		ctrl.header_freq,
		ctrl.report_precision,
		ctrl.maxit,
		ctrl.contract_factor,
		ctrl.expand_factor,
		ctrl.contract_threshold,
		ctrl.expand_threshold_rad,
		ctrl.expand_threshold_ap,
		ctrl.function_scale_factor,
		ctrl.precond_refresh_freq,
		ctrl.precond_ID,
//		ctrl.quasi_newton_method,
		ctrl.trust_iter
		);


	opt.run();

	const size_t Nparams = parameters.size();
	Tvec P(Nparams);
	Tvec grad(Nparams);
	Eigen::SparseMatrix<double> H(Nparams, Nparams);
	H.reserve(funObj.Htemplate.nonZeros());
// = funObj.Htemplate.cast<double>();

	double fval, radius;
	int iterations;
	MB_Status status;


	status = opt.get_current_state(P, fval, grad, H, iterations, radius);

	// get log determinant of hessian
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
	ldlt.compute(H);
	if (ldlt.info() != Eigen::Success) {
		Rcpp::warning("LLT factorization failed; H may not be SPD");
	}

	// log likelhood

	const auto& D = ldlt.vectorD();
	double logdetH = 0.0;
	for (int k = 0; k < D.size(); ++k) {
		logdetH += std::log(D[k]);
	}
	const double halfLogDet = logdetH/2;

	const double minusLogLik = fval + halfLogDet +
 	 Nparams * 0.918938533204672669541; // format(log(2*pi)/2, digits=22) 

 	 // save chol
	const auto& Pidx = ldlt.permutationP().indices();    
 	const Eigen::SparseMatrix<double> L = ldlt.matrixL();

	Rcpp::NumericVector D_R(D.size());
  	for (int k = 0; k < D.size(); ++k) D_R[k] = D[k];

  	Rcpp::IntegerVector P_R(Pidx.size());
  	for (int k = 0; k < Pidx.size(); ++k) P_R[k] = Pidx[k];  

 	Rcpp::List cholHessian = Rcpp::List::create(
		Rcpp::Named("P") = P_R,
		Rcpp::Named("D") = D_R,
		Rcpp::Named("L") = eigen_to_dgC(L, true)	// true for lower triangle only
 		);

	 Rcpp::NumericVector solution(P.size());
	for(size_t D=0;D<solution.size();D++) {
		solution[D] = P[D];
	}

 	 Rcpp::List res = Rcpp::List::create(
 	 	Rcpp::Named("minusLogLik") = Rcpp::wrap(minusLogLik),
 	 	Rcpp::Named("fval") = Rcpp::wrap(fval),
 	 	Rcpp::Named("halfLogDet") = Rcpp::wrap(halfLogDet),
 	 	Rcpp::Named("solution") = Rcpp::wrap(solution),
 	 	Rcpp::Named("gradient") = Rcpp::wrap(grad),	
 	 	Rcpp::Named("hessian") = eigen_to_dgC(H),
 	 	Rcpp::Named("cholHessian") = cholHessian,
 	 	Rcpp::Named("iterations") = Rcpp::wrap(iterations),
 	 	Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
 	 	Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
 	 	Rcpp::Named("method") = Rcpp::wrap("Sparse")
 	 	);

	if(!adPackFull) {
 	 return(res);		
	}

	if(config.verbose) {
		Rcpp::Rcout << "computing full hess";
	}

	const size_t Nbeta = config.beta.size(), Ntheta = config.theta.size();
	Eigen::VectorXd fullParameters(Nparams + Nbeta + Ntheta);
	for(size_t D=0;D<Nbeta;D++) {
		fullParameters[D] = config.beta[D];
	}
	for(size_t D=0;D<Nparams;D++) {
		fullParameters[D+Nbeta] = P[D];
	}
	for(size_t D=0;D<Ntheta;D++) {
		fullParameters[D+Nbeta+Nparams] = config.theta[D];
	}

	auto hessianIPLowerHere = config.hessianIPLower_outer;

	AD_Func_Opt funObjOuter(*adPackFull, hessianIPLowerHere);
	if(config.verbose) {
		Rcpp::Rcout << ".";
	}


	Eigen::SparseMatrix<double> hessianFull = funObj.Htemplate.cast<double>();
	Eigen::VectorXd gradFull(fullParameters);
	double f;

	if(config.verbose) {
		Rcpp::Rcout << ".";
	}



	funObjOuter.get_fdfh(fullParameters, f, gradFull, hessianFull);
	if(config.verbose) {
		Rcpp::Rcout << ".\n";
	}

 	 res["hessian_full"] = eigen_to_dgC(hessianFull);
	 res["gradient_full"] = Rcpp::wrap(gradFull);	

 	 return(res);
 	}

#endif
