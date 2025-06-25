#include"hpol.hpp"
#include<omp.h>


// [[Rcpp::export]]
Rcpp::List test3(
  Rcpp::NumericVector X,
  Rcpp::NumericVector U,
  Rcpp::NumericVector W){

  std::vector<double> w = Rcpp::as<std::vector<double>>(W);
  std::vector<double> x = Rcpp::as<std::vector<double>>(X);
  std::vector<double> direction = Rcpp::as<std::vector<double>>(U);
  int Nparams = X.size();

      Rcpp::Rcout << "Nparms" << Nparams << "\n";

  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = X[D];  // Initialize CppAD variables
  }
    CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y(1,0.0);
  for (size_t D = 0; D < Nparams; D++) {
    y[0] += ad_params[D]*ad_params[D]*ad_params[D];
  }   

  Rcpp::Rcout << "b\n";

  CppAD::ADFun<double> fun(ad_params, y);
        Rcpp::Rcout << "b1\n";

  std::vector<double> y_val(1);
  Rcpp::Rcout << "b2\n";

  y_val = fun.Forward(0, x);
          Rcpp::Rcout << "b3\n";
  fun.Forward(1, direction);
        Rcpp::Rcout << "c\n";


  auto taylor3 = fun.Reverse(3, w);

        Rcpp::Rcout << "d\n";

  Rcpp::NumericVector result(taylor3.size());
  for(size_t D =0; D<result.size();++D) result[D] = taylor3[D];
        Rcpp::Rcout << "e\n";

  return Rcpp::List::create(
    Rcpp::Named("taylor3") = result
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

  Rcpp::List sparsityR = config["sparsity"];
  Rcpp::IntegerVector HrowR = sparsityR["i"]; 
  Rcpp::IntegerVector HcolR = sparsityR["j"];
  Rcpp::IntegerVector HpR = sparsityR["p"];

  const Rcpp::S4 A(data["ATp"]), X(data["XTp"]);
  const Rcpp::IntegerVector dimsX = X.slot("Dim");
  const Rcpp::IntegerVector dimsA = A.slot("Dim");

  const int Nparams = parameters.size();
  const int Nbeta = dimsX[0];
  const int Ngamma = dimsA[0];
  const int Ntheta = Nparams - Nbeta - Ngamma;
  const int NbetaTheta = Nbeta + Ntheta;
  const int Nnonzero = HrowR.size();


  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);
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

 std::vector<double> w(3, 0.0); 
  w[2] = 1.0;



  // gammaParamMat = expand.grid(seq(Nbeta, Ngamma), c(seq(0, Nbeta-1), seq(Nbeta+Ngamma, len=Ntheta)))
 // combinations of gammas and thetas to compute
 const int NgammaCombo = Ngamma * NbetaTheta;
 Rcpp::IntegerMatrix gammaPramMat(NgammaCombo, 4);
 Rcpp::CharacterVector therownames = 
  Rcpp::CharacterVector::create(
    "gamma", "betaTheta", "gammaInParams", "betaThetaInParams"
  );
 gammaPramMat.attr("dimnames") = Rcpp::List::create(R_NilValue, therownames);

 int idx = 0;
 for (int Dgamma = 0; Dgamma < Ngamma; ++Dgamma) {
  for (int DbetaTheta = 0; DbetaTheta < NbetaTheta; ++DbetaTheta) {
    gammaPramMat(idx,0) = Dgamma;
    gammaPramMat(idx,1) = DbetaTheta;    
    gammaPramMat(idx,2) = Dgamma + Nbeta;// DgammaInParams
    gammaPramMat(idx,3) = DbetaTheta + (Ngamma) * (DbetaTheta >= Nbeta); //  DbetaThetaInParams 
    idx++;
  }
 }

Rcpp::NumericMatrix HgammaTheta(Ngamma, NbetaTheta);
Rcpp::NumericMatrix TgammaGammaTheta(HrowR.size(), NbetaTheta);
// a 3-D array Ngamma by NbetaTheta by NbetaTheta
Rcpp::NumericMatrix TgammaThetaTheta(Ngamma*NbetaTheta, NbetaTheta);

    if (verbose ) {
      Rcpp::Rcout << "starting 3rd deriv\n";
    }
    omp_set_num_threads(num_threads);

//  computing  T_iik + T_jjk + 2 T_ijk for all k.
// i's are gammas, j's beta or theta , k is either
  #pragma omp parallel
    {
  const int tid = omp_get_thread_num();

  #pragma omp for
 for(int DgammaCombo=0; DgammaCombo < 1;//NgammaCombo;
   ++DgammaCombo) {
    if (verbose ) {
      Rcpp::Rcout << "t " << tid << " i " << DgammaCombo << "; ";
    }



  int Dgamma = gammaPramMat(DgammaCombo,0);
  int DbetaTheta = gammaPramMat(DgammaCombo,1);
  int DgammaInParams = gammaPramMat(DgammaCombo,2);
  int DbetaThetaInParams = gammaPramMat(DgammaCombo,3);
  int Nnonzero = HpR[Dgamma+1] - HpR[Dgamma];

  std::vector<double> direction(Nparams, 0.0);
  direction[DgammaInParams]  = 1.0;     
  direction[DbetaThetaInParams]  = 1.0;     

  fun_threads[tid].Forward(1, direction);


  // first col is gradient
  auto taylor3 = fun_threads[tid].Reverse(3, w);

      if (verbose ) {
      Rcpp::Rcout << "size " << taylor3.size() << "\n";
    }

  // 2nd col H direction
  // saving H_gamma,gamma + H_theta,gamma
  HgammaTheta(Dgamma, DbetaTheta) = taylor3[Nparams + DgammaInParams];

  // 3rd col 0.5 * T_gamma,theta,k + T_gamma,gamma,k + T_theta,theta,k
  // Dgamma is column, take values of k which are gammas
  for (int Dnonzero = 0, DrowIndex = HpR[Dgamma]; Dnonzero < Nnonzero; ++Dnonzero, ++DrowIndex) {
    int Drow = HrowR[DrowIndex];
    TgammaGammaTheta(DrowIndex, DbetaTheta) = taylor3[2 * Nparams + Nbeta + Drow];
  }
  // now k's that are thetas or betas, put in TgammaThetaTheta[Dgamma, Dtheta2, Dtheta]
  for(int DbetaTheta2 = 0; DbetaTheta2 < NbetaTheta; DbetaTheta2++) {
    int DbetaTheta2InParams = DbetaTheta2 + (Ngamma) * (DbetaTheta2 >= Nbeta);
    TgammaThetaTheta(Dgamma + Ngamma * DbetaTheta2, DbetaTheta) = 
        taylor3[2 * Nparams + DbetaTheta2InParams];
  }
} // for
} // parallel


return Rcpp::List::create(
  Rcpp::Named("HgammaTheta") = HgammaTheta,
  Rcpp::Named("TgammaGammaTheta") = TgammaGammaTheta,
  Rcpp::Named("TgammaThetaTheta") = TgammaThetaTheta,
  Rcpp::Named("gammaPramMat")  = gammaPramMat
  );
}
