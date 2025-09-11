#include"hpol.hpp"
#include<omp.h>


#define DEBUG

#ifdef DEBUG

//' @export
// [[Rcpp::export]]
Rcpp::List test3(Rcpp::NumericVector X, Rcpp::NumericVector U, Rcpp::NumericVector V, Rcpp::NumericVector W){

  std::vector<double> w = Rcpp::as<std::vector<double>>(W);
  std::vector<double> x = Rcpp::as<std::vector<double>>(X);
  std::vector<double> direction1 = Rcpp::as<std::vector<double>>(U);
  std::vector<double> direction2 = Rcpp::as<std::vector<double>>(V);
  int Nparams = X.size();

  Rcpp::Rcout << "Nparms" << Nparams << "\n";

  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = X[D];  // Initialize CppAD variables
  }
    CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

    CppAD::vector<CppAD::AD<double>> y(1,0.0);
  for (size_t D = 0; D < Nparams; //Nparams; 
    D++) {
    y[0] += ad_params[D]*(1+D);
}   
y[0] = y[0]*y[0]*y[0]*y[0];

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
#endif

/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::List derivForLaplace(
  Rcpp::NumericVector parameters, // beta, gamma, theta
  Rcpp::List data, 
  Rcpp::List config // sparsity$second full and random sparsity$third jk and ijk
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
 const int Ngamma = dimsA[0];
 const int Ntheta = Nparams - Nbeta - Ngamma;

// hessian of random effects, lower triangle only, column format
 Rcpp::List sparsity = config["sparsity"];
 Rcpp::List hessianList = sparsity["second"];
 Rcpp::List hessian = hessianList["full"];
 Rcpp::IntegerVector Hrow = hessian["i"]; 
 Rcpp::IntegerVector Hp = hessian["p"];

// sparsity for third deriv
 Rcpp::List thirdList = sparsity["third"];
 Rcpp::List sparsityThirdIjk = thirdList["ijk"];   //  "p" "i" "j" "k"
 Rcpp::List sparsityThirdJk = thirdList["jk"];    // "pEnd" "n"    "p"    "j"    "k" 

 Rcpp::IntegerVector sparsityJkJ = sparsityThirdJk["j"];
 Rcpp::IntegerVector sparsityJkK = sparsityThirdJk["k"];
 Rcpp::IntegerVector sparsityJkP = sparsityThirdJk["p"];
 Rcpp::IntegerVector sparsityJkPend = sparsityThirdJk["pEnd"];
 Rcpp::IntegerVector sparsityIjkI = sparsityThirdIjk["i"];

 const int Nthird = sparsityIjkI.size();
 const int Npairs = sparsityJkJ.size();

// output
 Rcpp::NumericVector hessianOut(Hrow.size());
 Rcpp::NumericVector diagOut(Hrow.size());
 Rcpp::NumericVector Tijk(Nthird);    


#ifdef DEBUG
 Rcpp::NumericMatrix hessianDense(Nparams, Nparams);
 Rcpp::NumericMatrix thirdDiagDense(Nparams, Nparams);
 Rcpp::NumericMatrix TijkDense(Nparams, Npairs);
#endif

// set up autodiff function
 CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
 for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = 
  objectiveFunctionInternal(ad_params, data, config);  
  CppAD::ADFun<double> fun(ad_params, y);

  std::vector<double> x_val(Nparams);
  std::vector<double> y_val(1);
  std::vector<double> w{0.0, 0.0, 1.0};  

  for (size_t i = 0; i < Nparams; ++i) {
    x_val[i] = parameters[i];
  }
  y_val = fun.Forward(0, x_val);

  // Replicate fun object for each thread
  std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
  for (int i = 0; i < num_threads; ++i) {
    fun_threads[i] = fun;
  }
  if (verbose ) {
    Rcpp::Rcout << "starting 3rd deriv, " << num_threads << " threads\n";
  }
  omp_set_num_threads(num_threads);

//  computing  T_iik + T_jjk + 2 T_ijk for all k.
// i's are gammas, j's beta or theta , k is either
  #pragma omp parallel
  {
    const int tid=omp_get_thread_num();

//    std::vector<double> direction(Nparams, 0.0);


  // diagonal entries to subtract off
  // computing T_kii and H_ki, columnns are i, rows are k
    #pragma omp for
    for(int i=0; i < Nparams; ++i) {

      const int HessianStart = Hp[i+1];
      const int HessianEnd = Hp[i+1];


//    std::fill(direction.begin(), direction.end(), 0.0);
      std::vector<double> direction(Nparams, 0.0);
      direction[i]  = 1.0;     


      fun_threads[tid].Forward(0, x_val);
      fun_threads[tid].Forward(1, direction);
      fun_threads[tid].Forward(2, direction);

      auto taylor3 = fun_threads[tid].Reverse(3, w);

  // fill in the T_iik and Hessian
      for(
        int Dj=HessianStart; 
        Dj<HessianEnd; 
        Dj++){
        int indexHere = 3*Hrow[Dj];
      // second column is hessian
      hessianOut[Dj] = taylor3[1+indexHere];
      // first column is T_ijj/2 + H_ij
      diagOut[Dj] = 2*(taylor3[indexHere] - hessianOut[Dj]);
    }



#ifdef DEBUG
    // store dense
    for(int Dj=0; Dj<Nparams; Dj++){
      hessianDense(Dj, i) = taylor3[1+3*Dj];
      thirdDiagDense(Dj, i) = taylor3[3*Dj];
    }
#endif  

 } // for i diagonal bit

// off diag T_ijk

    #pragma omp for
 for(int Djk=0; Djk < Npairs; ++Djk) {
  const int Dj = sparsityJkJ[Djk];
  const int Dk = sparsityJkK[Djk];

  std::vector<double> direction(Nparams, 0.0);
  direction[Dk]  = direction[Dj] = 1.0;     

  fun_threads[tid].Forward(0, x_val);
  fun_threads[tid].Forward(1, direction);
  fun_threads[tid].Forward(2, direction);

  // first column is third deriv combination
  //  T_iik + T_jjk + 2 T_ijk 
  auto taylor3 = fun_threads[tid].Reverse(3, w);


  for(int Di=sparsityJkP[Djk]; Di<sparsityJkPend[Djk]; Di++){
    int indexHere = 3*sparsityIjkI[Di];
    Tijk[Di] = taylor3[indexHere]; // /to do subttract off Tiik, Tjjk
  }

#ifdef DEBUG
  for(int Di=0; Di<Nparams; Di++){
    TijkDense(Di, Djk) = taylor3[3*Di];
  }
#endif  

 } // Djk


} // parallel


Rcpp::List resultList = Rcpp::List::create(
  Rcpp::Named("third") = Tijk,
  Rcpp::Named("diag") = diagOut,
  Rcpp::Named("hessian") = hessianOut
#ifdef DEBUG
  ,
  Rcpp::Named("denseHessian") = hessianDense,
  Rcpp::Named("denseDiag") = thirdDiagDense,
  Rcpp::Named("denseThird") = TijkDense
#endif  
  );


return resultList;
}
