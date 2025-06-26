#include"hpol.hpp"
#include<omp.h>


// [[Rcpp::export]]
Rcpp::List test3(
  Rcpp::NumericVector X,
  Rcpp::NumericVector U,
  Rcpp::NumericVector V,
  Rcpp::NumericVector W){

  std::vector<double> w = Rcpp::as<std::vector<double>>(W);
  std::vector<double> x = Rcpp::as<std::vector<double>>(X);
  std::vector<double> direction1 = Rcpp::as<std::vector<double>>(U);
  std::vector<double> direction2 = Rcpp::as<std::vector<double>>(V);
  int Nparams = X.size();

      Rcpp::Rcout << "Nparms" << Nparams << "\n";

  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < 2; D++) {
    ad_params[D] = X[D];  // Initialize CppAD variables
  }
    CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y(1,0.0);
  for (size_t D = 0; D < 2; //Nparams; 
    D++) {
    y[0] += ad_params[D]*ad_params[D]*ad_params[D];
  }   

  Rcpp::Rcout << "b\n";

  CppAD::ADFun<double> fun(ad_params, y);
        Rcpp::Rcout << "b1\n";

  std::vector<double> y_val(1);
  Rcpp::Rcout << "b2\n";

  y_val = fun.Forward(0, x);
          Rcpp::Rcout << "b3\n";
  auto taylor1 = fun.Forward(1, direction1);
  auto taylor2 = fun.Forward(2, direction2);
        Rcpp::Rcout << "c\n";


  auto taylor3 = fun.Reverse(3, w);

        Rcpp::Rcout << "d\n";

  Rcpp::NumericVector result(taylor3.size()), result2(taylor2.size());
  for(size_t D=0;D<result.size();D++) result[D] = taylor3[D];
  for(size_t  D=0;D<result2.size();D++) result2[D] = taylor2[D];

  return Rcpp::List::create(
    Rcpp::Named("taylor3") = result,
    Rcpp::Named("taylor2") = result2
  );

}

// [[Rcpp::export]]
Rcpp::List derivForLaplace(
  Rcpp::NumericVector parameters, // gamma, beta, theta
  Rcpp::List data, 
  Rcpp::List config
  ) {

	int num_threads = 1;
	if (config.containsElementNamed("num_threads"))
		num_threads = Rcpp::as<int>(config["num_threads"]);

  bool verbose = false;
   if (config.containsElementNamed("verbose")) {
   verbose = Rcpp::as<bool>(config["verbose"]);
   }


  const Rcpp::S4 A(data["ATp"]), X(data["XTp"]);
  const Rcpp::IntegerVector dimsX = X.slot("Dim");
  const Rcpp::IntegerVector dimsA = A.slot("Dim");

  const int Nparams = parameters.size();
  const int Nbeta = dimsX[0];

// hessian of random effects, lower triangle only, column format
  Rcpp::List hessianRandom = config["sparsity"];
  Rcpp::IntegerVector Hrow = hessianRandom["i"]; 
  Rcpp::IntegerVector Hp = hessianRandom["p"];

// indices for third deriv tensor 
  Rcpp::DataFrame parametersGamma = 
    Rcpp::as<Rcpp::DataFrame>(config["parametersGamma"]);
  Rcpp::IntegerVector PinThird = parametersGamma["PinThird"];
  Rcpp::IntegerVector NinHessian = parametersGamma["NinHessian"];
  Rcpp::IntegerVector PinHessian = parametersGamma["PinHessian"];
  Rcpp::IntegerVector gammaInParams = parametersGamma["gammaInParams"];
  Rcpp::IntegerVector paramInParams = parametersGamma["paramInParams"];

  int NgammaCombo = parametersGamma.nrow();
// to store the third deriv tensor
  const int Nthird = 
    PinThird[NgammaCombo-1] +
      NinHessian[NgammaCombo-1];
  Rcpp::NumericVector result(Nthird);    


// indices computation of entries T_iij
  Rcpp::List indexForDiag = config["indexForDiag"];
  Rcpp::IntegerVector indexForDiagI = indexForDiag["i"];
  Rcpp::IntegerVector indexForDiagP = indexForDiag["p"];
  const int Ndiag = indexForDiagP.size()-1;
  Rcpp::NumericVector forDiag(indexForDiagI.size());



  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = 
    objectiveFunctionInternal(ad_params, data, config);
  CppAD::ADFun<double> fun(ad_params, y);
  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) {
    x_val[i] = parameters[i];
  }

  std::vector<double> y_val(1);
  y_val = fun.Forward(0, x_val);

    // Replicate fun object for each thread
  std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
  for (int i = 0; i < num_threads; ++i) fun_threads[i] = fun;

  std::vector<double> w{0.0, 0.0, 1.0};  


    if (verbose ) {
      Rcpp::Rcout << "starting 3rd deriv, " << num_threads << " threads\n";
    }
    omp_set_num_threads(num_threads);

//  computing  T_iik + T_jjk + 2 T_ijk for all k.
// i's are gammas, j's beta or theta , k is either
  #pragma omp parallel
    {
  const int tid=omp_get_thread_num();

  #pragma omp for
 for(int DgammaCombo=0; DgammaCombo < NgammaCombo;   ++DgammaCombo) {

  const int DgammaInParams = gammaInParams[DgammaCombo];
  const int DbetaThetaInParams = paramInParams[DgammaCombo];
  const int Nnonzero = NinHessian[DgammaCombo];

  std::vector<double> direction(Nparams, 0.0);
  direction[DgammaInParams]  = 1.0;     
  direction[DbetaThetaInParams]  = 1.0;     

  fun_threads[tid].Forward(0, x_val);
  fun_threads[tid].Forward(1, direction);
  fun_threads[tid].Forward(2, direction);


  // result is row-wise, first column is third deriv combination
  auto taylor3 = fun_threads[tid].Reverse(3, w);

  const int PinHessianHere = PinHessian[DgammaCombo];
  const int PinThirdHere = PinThird[DgammaCombo];
  for(int Dthird=0; Dthird < Nnonzero; Dthird++)  {
    int Drow = Hrow[PinHessianHere + Dthird] + Nbeta;
    result[PinThirdHere+Dthird] = taylor3[3*Drow];
  }
 } // for
} // parallel

    if (verbose ) {
      Rcpp::Rcout << "diag parts, " << num_threads << " threads\n";
    }

  #pragma omp parallel
    {
  const int tid=omp_get_thread_num();

  // diagonal entries to subtract off
  // computing T_iik,  columnns are i, rows are k
    #pragma omp for
  for(int i=0; i < Ndiag; ++i) {
    const int Pstart = indexForDiagP[i];
    const int Nhere = indexForDiagP[i+1] - Pstart;

    if(Nhere == 0) continue;

    std::vector<double> direction(Nparams, 0.0);
    direction[i]  = 1.0;     

  fun_threads[tid].Forward(0, x_val);
  fun_threads[tid].Forward(1, direction);
  fun_threads[tid].Forward(2, direction);

  auto taylor3 = fun_threads[tid].Reverse(3, w);

  for(int Dj=0; Dj<Nhere; Dj++){
      const int indexHere = Pstart+Dj;
      forDiag[indexHere] = 
        taylor3[3*indexForDiagI[indexHere]];
  }

 } // for
} // parallel


return Rcpp::List::create(
  Rcpp::Named("result") = result,
  Rcpp::Named("forDiag") = forDiag
  );
}
