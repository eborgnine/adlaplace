// funcitons useful for debugging in R.  not needed for model fitting

#ifndef DEBUGGING_HPP
#define DEBUGGING_HPP


double jointLogDens_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	Config configC(config);
	double result = 0.0;

	if (adPack == R_NilValue) {
		return(result);
	}
    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}
	const bool inner = (Nparams == configC.start_gamma.size());
	auto hessianIPLowerHere = inner?configC.hessianIPLower_inner:configC.hessianIPLower_outer;


	AD_Func_Opt funObj(*xp, hessianIPLowerHere);

	Eigen::VectorXd parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	cppad_parallel_setup(configC.num_threads);

	if(configC.verbose) {
		Rcpp::Rcout << "d";
	}


	funObj.get_f(parametersC, result);

	if(configC.verbose) {
		Rcpp::Rcout << "e " << result << "\n";
	}


	return(result);

}


Rcpp::NumericVector grad_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{
	Config configC(config);

	if (adPack == R_NilValue) {
		return(Rcpp::NumericVector());
	}

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}
	const bool inner = (Nparams == configC.start_gamma.size());
	auto hessianIPLowerHere = inner?configC.hessianIPLower_inner:configC.hessianIPLower_outer;

	AD_Func_Opt funObj(*xp, hessianIPLowerHere);

	Rcpp::NumericVector result(Nparams);


	Eigen::VectorXd parametersC(Nparams);
	Eigen::VectorXd resultC(Nparams);

	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	cppad_parallel_setup(configC.num_threads);

	funObj.get_df(parametersC, resultC);

	for(size_t D=0; D<Nparams;D++) {
		result[D] = resultC[D];
	}
	return(result);

}

Rcpp::S4 hessian_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	Rcpp::List config)
{

	if (adPack == R_NilValue) {
		return(R_NilValue);
	}

	Config configC(config);

    // Rehydrate external pointer
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}

	const bool inner = (Nparams == configC.start_gamma.size());
	auto hessianIPLowerHere = inner?configC.hessianIPLower_inner:configC.hessianIPLower_outer;

	AD_Func_Opt funObj(*xp, hessianIPLowerHere);


	Eigen::VectorXd parametersC(Nparams);
	
	for(size_t D=0; D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	Eigen::SparseMatrix<double> resultC = funObj.Htemplate.cast<double>();

	cppad_parallel_setup(configC.num_threads);

	if(configC.verbose) {
		Rcpp::Rcout << "hessian params " << parametersC.size() << " groups "	<< (*xp).size();
	}

	funObj.get_hessian(parametersC, resultC);
	if(configC.verbose) {
		Rcpp::Rcout  << "\n";
	}

	resultC.makeCompressed();

/*	Eigen::SparseMatrix<double> Hfull(
		resultC.selfadjointView<Eigen::Lower>()
		);*/
	Rcpp::S4 out = eigen_to_dgC(resultC);

	return(out);
}


double jointLogDensNoAdfun_backend(
	Rcpp::NumericVector parameters, 
	Rcpp::List dataR,
	Rcpp::List configR)
{

	Config config(configR);
	Data data(dataR);


	const size_t Nparams = parameters.size();
	const size_t Nbeta = config.beta.size();
//	const size_t Ntheta = config.theta.size();
	const size_t Ngamma = config.start_gamma.size();

	size_t Ngroups = config.groups.ncol();
	if(Ngroups==0) {
		// groups not provided, check for elgm matrix
		Ngroups = data.elgm_matrix.ncol();
		if(Ngroups == 0) {
			// no elgm matrix, use y
			Ngroups = data.y.size();
		}
	}

	const bool inner = (Nparams == Ngamma);

	CppAD::vector<CppAD::AD<double>> parametersC(Nparams);
	for(size_t D=0; D<Nparams;D++) {
		double parHere = parameters[D];
		parametersC[D] = parHere;
	}

	CppAD::AD<double> dataPart=0.0, extraPart, randomPart;
	if(inner) {
		for(size_t D=0;D<Ngroups;D++) {
			if(config.verbose) {
				Rcpp::Rcout << D << " ";
			}
			auto dataPartHere= logDensObs(parametersC, config.beta, config.theta, data, config, D);
			dataPart += dataPartHere[0];
		}
		auto randomHere = logDensRandom(
			parametersC,
			config.theta, data, config);   
		randomPart = randomHere[0];

		auto extraHere = logDensExtra(
			config.theta, data, config);   
		extraPart = extraHere[0];
	} else {
		for(size_t D=0;D<Ngroups;D++) {
			auto dataPartHere= logDensObs(
				slice(parametersC, Nbeta, Nbeta + Ngamma), //gamma
								slice(parametersC, 0, Nbeta), // beta
								slice(parametersC, Nbeta + Ngamma, Nparams), //theta
								data, config, D);
			dataPart += dataPartHere[0];
		}
		auto randomHere = logDensRandom(
				slice(parametersC, Nbeta, Nbeta + Ngamma), //gamma
												slice(parametersC, Nbeta + Ngamma, Nparams), //theta
												data, config);   
		randomPart = randomHere[0];

		auto extraHere = logDensExtra(
			slice(parametersC, Nbeta + Ngamma, Nparams), //theta
			data, config);   
		extraPart = extraHere[0];

	}
	CppAD::AD<double> result1 = dataPart + randomPart + extraPart;



	double result = CppAD::Value(result1);

	if(config.verbose) {
		Rcpp::Rcout << "Ngroups " << Ngroups << "inner " << inner <<  " d " << dataPart << " r " << randomPart << " e " 
			<< extraPart << " result1 " << result1 << " result " << result << "\n";
	}


	return(result);

}




#endif
