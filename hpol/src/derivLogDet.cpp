

Rcpp::List derivForLaplace(
  Rcpp::NumericVector parameters, // gamma, beta, theta
  Rcpp::List data, 
  Rcpp::List config
  ) {

	int num_threads = 1;
	if (config.containsElementNamed("num_threads"))
		num_threads = Rcpp::as<int>(config["num_threads"]);

	const Rcpp::S4 A(data["ATp"]), X(data["XTp"]);
	const Rcpp::IntegerVector dimsX = X.slot("Dim"), dimsA = A.slot("Dim"),
	const int Nparams = parameters.size(),
		Nbeta =dimsX[0], Ngamma = dimsA[0], 
		Ntheta = Nparams - Nbeta-Ngamma, NbetaTheta = Nbeta+Ntheta;

	Rcpp::List sparsityR = config["sparsity"];
	Rcpp::IntegerVector HrowR = sparsityR["i"]; 
	Rcpp::IntegerVector HpR = sparsityR["p"];
	Rcpp::NumericMatrix HgammaParam(Ngamma, NbetaTheta);
	Rcpp::NumericMatrix thirdDeriv(HrowR.size(), NbetaTheta);

	CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
	for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
}
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);
  CppAD::ADFun<double> fun(ad_params, y);
  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) x_val[i] = parameters[i];

  	std::vector<double> y_val(1);
  y_val = fun.Forward(0, x_val);

    // Replicate fun object for each thread
  std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
  for (int i = 0; i < num_threads; ++i) fun_threads[i] = fun;


  	#pragma omp parallel
  	for(int Dgamma=0; Dgamma < Ngamma; ++Dgamma) {
  		int pHere = HpR[Dgamma];   
  		int NpHere = HpR[Dgamma+1] - pHere;
  		const int tid = thread_number();

    // Create a direction vector for Forward(1)
  		std::vector<double> direction(Nparams, 0.0);
    direction[Dgamma] = 1.0; // derivative with respect to parameter Dparam


    // Compute Taylor coefficients for first derivative
    auto taylor1 = fun_threads[tid].Forward(1, direction);


    for(int DbetaTheta=0;DbetaTheta<NbetaTheta;++DbetaTheta) {

    	std::vector<double> direction2(Nparams, 0.0);
    	direction2[Dgamma] = 1.0;
    	direction2[Dgamma + DbetaTheta] = 1.0;

    	auto taylor2 = fun_threads[tid].Forward(2, direction2);

    	HgammaParam(Dgamma, DbetaTheta) = taylor2;

    	for(int Dinner=0; Dinner < NpHere; ++Dinner) {
    		int Dp = pHere + Dinner;

    		std::vector<double> direction3(Nparams, 0.0);
    		direction3[Dgamma] = 1.0;
    		direction3[HrowR[Dp]] = 1.0;
    		direction3[Dgamma + DbetaTheta] = 1.0;
    		auto taylor3 = fun_threads[tid].Forward(3, direction3);

    		thirdDeriv(Phere+Dinner, DbetaTheta) = taylor3;
    	}

    }

}