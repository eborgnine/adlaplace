#ifndef ADFUN_HPP
#define ADFUN_HPP


template <class Type>
CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &beta,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	);

template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	);

template <class Type>
CppAD::vector<CppAD::AD<double>> logDensRandom(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	);

inline std::vector<GroupPack> getAdFunOuter(
	const Data& data,
	const Config& config) {

	size_t Ngroups = config.groups.ncol();
	if(Ngroups==0) {
		// groups not provided, check for elgm matrix
		Ngroups = data.elgm_matrix.ncol();
		if(Ngroups == 0) {
			// no elgm matrix, use y
			Ngroups = data.y.size();
		}
	}


	std::vector<GroupPack> result(Ngroups+2L);

		const size_t Ngamma = config.start_gamma.size();
		const size_t Nbeta = config.beta.size();
		const size_t Ntheta = config.theta.size();
		const size_t Nparams = Nbeta + Ngamma + Ntheta;

	if(config.verbose) {
		Rcpp::Rcout << "outer, groups " << Ngroups << " Nbeta " << Nbeta << " Ntheta " <<
		 Ntheta << " Ngamma " << Ngamma << " Nparams " << Nparams << "\n";
	}

		CppAD::vector<CppAD::AD<double>> ad_params_G(Nparams);
		for(size_t D=0;D<Nbeta;D++) {
			ad_params_G[D] = config.beta[D];
		}
		for(size_t Dgamma=0;Dgamma<Ngamma;Dgamma++) {
			ad_params_G[Dgamma+Nbeta] = config.start_gamma[Dgamma];
		}
		for(size_t Dtheta=0;Dtheta<Ntheta;Dtheta++) {
				ad_params_G[Dtheta+Nbeta+Ngamma] = config.theta[Dtheta];
		} 

	if(config.verbose) {
		Rcpp::Rcout << ".";
	}
# pragma omp parallel
	{
			CppAD::vector<CppAD::AD<double>> ad_params = ad_params_G;

    # pragma omp for schedule(dynamic,1) nowait
		for(size_t D=0;D<Ngroups;D++) {


			CppAD::Independent(ad_params);
			auto resultHere = logDensObs(
				slice(ad_params, Nbeta, Nbeta + Ngamma), // gamma
				slice(ad_params, 0, Nbeta), // beta
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config, D);
			CppAD::ADFun<double> fun(ad_params, resultHere);

			result[D].fun = std::move(fun);
		}


# pragma omp single nowait
		{
			CppAD::Independent(ad_params);

			auto resultHere = logDensRandom(
				slice(ad_params, Nbeta, Nbeta + Ngamma), // gamma
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[Ngroups].fun = std::move(fun);

		}

		# pragma omp single nowait
		{
			CppAD::Independent(ad_params);

			auto resultHere = logDensExtra<CppAD::AD<double>>(
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[Ngroups+1].fun = std::move(fun);

		}

	} // parallel

	if(config.group_sparsity_outer.size()) {
		if(config.verbose) {
			Rcpp::Rcout << "add sparse\n";
		}
		addSparsityToAdFun(result, config.group_sparsity_outer, config.verbose);
	} else {
		if(config.verbose) {
			Rcpp::Rcout << "no sparse pattern provided\n";
		}
	}

	return(result);
}


inline std::vector<GroupPack> getAdFunInner(
	const Data& data,
	const Config& config) {

	auto parameters = config.start_gamma;

	size_t Ngroups = config.groups.ncol();
	if(Ngroups==0) {
		// groups not provided, check for elgm matrix
		Ngroups = data.elgm_matrix.ncol();
		if(Ngroups == 0) {
			// no elgm matrix, use y
			Ngroups = data.y.size();
		}
	}

	std::vector<GroupPack> result(Ngroups+2L);

	if(config.verbose) {
		Rcpp::Rcout << " adfun Nparams " << parameters.size() << " groups " << Ngroups << " ";
	}

# pragma omp parallel
	{
		const size_t Nparams = parameters.size();
		CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
		for(size_t D=0;D<Nparams;D++) {
			ad_params[D] = parameters[D];
		}


    # pragma omp for schedule(dynamic,1) nowait
		for(size_t D=0;D<Ngroups;D++) {

			CppAD::Independent(ad_params);   
			auto resultHere = logDensObs(
				ad_params, 
				config.beta, 
				config.theta,
				data, config, D);
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);
		}



# pragma omp single nowait
		{
			CppAD::Independent(ad_params);
			auto resultHere = logDensRandom(ad_params, 
				config.theta, data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[Ngroups].fun = std::move(fun);
		}

		# pragma omp single nowait
		{
			// adds the extra part to the likelihood
			// won't contribute to derivatives
			CppAD::Independent(ad_params);
			CppAD::vector<double> resultExtra = logDensExtra(
				config.theta, 
				data, config);   
			CppAD::vector<CppAD::AD<double>> resultHere(1, 0.0);
			resultHere[0] += resultExtra[0];
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[Ngroups+1].fun = std::move(fun);

		}

	} // parallel

		if(config.verbose) {
			Rcpp::Rcout << "done adfun inner\n";
		}


	if(config.group_sparsity_inner.size()) {
		if(config.verbose) {
			Rcpp::Rcout << "add sparse " <<  config.group_sparsity_inner.size() << "\n";
		}
		addSparsityToAdFun(result, config.group_sparsity_inner);
	} else {
		if(config.verbose) {
			Rcpp::Rcout << "no sparse pattern provided\n";
		}
	}

	return(result);
}

inline double jointLogDensNoAdfun_backend(
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
