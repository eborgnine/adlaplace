#ifndef RINTERFACES_BACKEND_HPP
#define RINTERFACES_BACKEND_HPP


SEXP getAdFun_backend(
	Rcpp::List data, 
	Rcpp::List config,
	const bool inner=false)
{

	Data dataC(data);
	Config configC(config);
	if(configC.verbose) {
		Rcpp::Rcout << "getting ad function, inner " << inner << " ";
	}

	cppad_parallel_setup(configC.num_threads);

	std::vector<GroupPack> adPack = inner?
	getAdFunInner(dataC, configC):getAdFunOuter(dataC, configC);
	auto* ptr = new std::vector<GroupPack>(std::move(adPack));
  Rcpp::XPtr<std::vector<GroupPack>> xp(ptr, /*deleteOnFinalizer=*/true);

	xp.attr("class") = "adpack_ptr";
	return xp;
}


Rcpp::NumericVector traceHinvT_backend(
	const Rcpp::NumericVector parameters,
	const Rcpp::S4& LinvPt,
	const Rcpp::S4& LinvPtColumns,
	const Rcpp::List config,
	SEXP adPack
	) {

	const Config configC(config);
	const size_t Nparams = parameters.size();
	CppAD::vector<double> parametersC(Nparams);


	for(size_t D=0; D<Nparams; D++) {
		parametersC[D] = parameters[D];
	}

	CSCMatrix LinvPtC = makeCSC(LinvPt);
	CSCMatrix LinvPtColumnsC = makeCSC(LinvPtColumns);

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);


	cppad_parallel_setup(configC.num_threads);

	if(configC.verbose) {
		Rcpp::Rcout << ".";
	}

	auto resultC = traceHinvT(parametersC, LinvPtC, LinvPtColumnsC, configC, *xp);
	if(configC.verbose) {
		Rcpp::Rcout << ".";
	}

	Rcpp::NumericVector result(Nparams);
	for(size_t D=0; D<Nparams; D++) {
		result[D] = resultC[D];
	}
	return(result);
}

Rcpp::List inner_opt_backend(
	Rcpp::NumericVector parameters, 
	Rcpp::List data,
	Rcpp::List control, 
	Rcpp::List config,
	SEXP adPackFull = R_NilValue)
{

	Config configC(config);
	Data dataC(data);
	TrustControl ctrl(control); 

    // Optional full pack pointer
	std::vector<GroupPack>* adPackFullPtr = nullptr;

	if (adPackFull != R_NilValue) {

		if(configC.verbose) {
			Rcpp::Rcout << "outer adpack provided";
		}

        // Construct a *temporary* XPtr from adPackFull…
		Rcpp::XPtr<std::vector<GroupPack>> xpFull(adPackFull);
        // …and grab the raw pointer. The memory is owned by R, not by xpFull,
        // so this pointer remains valid after xpFull goes out of scope.
		adPackFullPtr = xpFull.get();
	} else {
		if(configC.verbose) {
			Rcpp::Rcout << "no outer adpack";
		}
	}

	using Tvec  = Eigen::VectorXd;

	const size_t Nparams = parameters.size();
	Tvec parametersC(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	cppad_parallel_setup(configC.num_threads);

	std::vector<GroupPack> adPack = getAdFunInner(dataC, configC);

	if(configC.verbose) {
		Rcpp::Rcout << "got inner adpack";
	}

	Rcpp::List res = inner_opt(parametersC, adPack, ctrl, configC, adPackFullPtr);

	return(res);
}


Rcpp::List inner_opt_adpack_backend(
	Rcpp::NumericVector parameters, 
	SEXP adPack,
	const Rcpp::List control, 
	const Rcpp::List config,
	SEXP adPackFull = R_NilValue)
{

	const Config configC(config);
	const TrustControl ctrl(control); 

	using Tvec  = Eigen::VectorXd;

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
    // Optional full pack pointer
	std::vector<GroupPack>* adPackFullPtr = nullptr;
	if (adPackFull != R_NilValue) {
        // Construct a *temporary* XPtr from adPackFull…
		Rcpp::XPtr<std::vector<GroupPack>> xpFull(adPackFull);
        // …and grab the raw pointer. The memory is owned by R, not by xpFull,
        // so this pointer remains valid after xpFull goes out of scope.
		adPackFullPtr = xpFull.get();
	}


	const size_t Nparams = parameters.size();
	const size_t NparamsFun = (*xp)[0].fun.Domain();
	if(Nparams != NparamsFun) {
		Rcpp::Rcout << "Warning parameters in " << Nparams << " parameters in AD fun " << NparamsFun << "\n";
	}

	Tvec parametersC(Nparams);
	for(size_t D=0;D<Nparams;D++) {
		parametersC[D] = parameters[D];
	}

	cppad_parallel_setup(configC.num_threads);

	Rcpp::List res = inner_opt(parametersC, *xp, ctrl, configC, adPackFullPtr);

	return(res);
}

Rcpp::List sparsity_backend(
	SEXP adPack,
	const Rcpp::NumericVector parameters,
	const bool verbose = false
	) {
	
	const size_t Nparams = parameters.size();
	CppAD::vector<double> parametersC(Nparams);
	for (size_t j = 0; j < Nparams; ++j)
		parametersC[j] = parameters[j];

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);

	Rcpp::List result = sparsity(*xp, parametersC, verbose);

	return(result);
}

Rcpp::List sparsity_backend(
	Rcpp::List data,
	Rcpp::List config
	) {

	Data dataC(data);
	Config configC(config);

	set_num_threads_wrapper(configC.num_threads);
	CppAD::thread_alloc::parallel_setup(
		configC.num_threads,
		[](){ return in_parallel_wrapper(); },
		[](){ return static_cast<size_t>(thread_num_wrapper()); }
		);

	auto adFun = getAdFunOuter(dataC, configC);

	if(configC.verbose) {
		Rcpp::Rcout << "..";
	}

	const size_t Ngamma = configC.start_gamma.size();
	const size_t Nbeta = configC.beta.size();
	const size_t Ntheta = configC.theta.size();
	const size_t Nparams = Nbeta + Ngamma + Ntheta;
	
	if(configC.verbose) {
		Rcpp::Rcout << ",";
	}
	CppAD::vector<double> ad_params(Nparams);
	for(size_t D=0;D<Nbeta;D++) {
		ad_params[D] = configC.beta[D];
	}
	for(size_t Dgamma=0;Dgamma<Ngamma;Dgamma++) {
		ad_params[Dgamma+Nbeta] = configC.start_gamma[Dgamma];
	}
	for(size_t Dtheta=0;Dtheta<Ntheta;Dtheta++) {
		ad_params[Dtheta+Nbeta+Ngamma] = configC.theta[Dtheta];
	} 
	if(configC.verbose) {
		Rcpp::Rcout << "getting sparse\n";
	}
	Rcpp::List result = sparsity(adFun, ad_params, configC.verbose);

	return(result);
}


#endif

