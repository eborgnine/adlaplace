#include"hpol.hpp"
#include<omp.h>


#define DEBUG



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



/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::List thirdDiagonals(
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
 bool dense = false;
 if (config.containsElementNamed("dense")) {
   dense = Rcpp::as<bool>(config["dense"]);
 }


 const int Nparams = parameters.size();

// output
// Rcpp::NumericVector gradient(Nparams);
 Rcpp::NumericMatrix hessianOut, thirdDiagOut, gradientOut;
 Rcpp::List sparsity, hessianList, hessian, diagSparsity;
 Rcpp::IntegerVector Hrow, Hp,  DiagRow, DiagP;

 if(dense){
  hessianOut = Rcpp::NumericMatrix(Nparams, Nparams);
  thirdDiagOut = Rcpp::NumericMatrix(Nparams, Nparams);
  gradientOut = Rcpp::NumericMatrix(Nparams, Nparams);
} else {
// hessian of random effects, lower triangle only, column format
  sparsity = config["sparsity"];
  hessianList = sparsity["second"];
  hessian = hessianList["full"];
  diagSparsity = hessianList["nonSymmetric"];

  Hrow = hessian["i"]; 
  Hp = hessian["p"];
  DiagRow = diagSparsity["i"]; 
  DiagP = diagSparsity["p"];
  
  hessianOut = Rcpp::NumericMatrix(DiagRow.size(), 1);
  thirdDiagOut = Rcpp::NumericMatrix(DiagRow.size(), 1);
  gradientOut  = Rcpp::NumericMatrix(Nparams, 1);
}

if (verbose ) Rcpp::Rcout << "objects allocated\n";

// set up autodiff function
CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);  
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

  omp_set_num_threads(num_threads);

  #pragma omp parallel
  {
    const int tid=omp_get_thread_num();
    fun_threads[tid].Forward(0, x_val);
    std::vector<double> direction(Nparams, 0.0), directionZeros(Nparams, 0.0);

  // diagonal entries to subtract off
  // computing T_kii and H_ki, columnns are i, rows are k
    #pragma omp for
    for(int Dk=0; Dk < Nparams; ++Dk) {

      std::fill(direction.begin(), direction.end(), 0.0);
      direction[Dk]  = 1.0;     

      fun_threads[tid].Forward(1, direction);
      fun_threads[tid].Forward(2, directionZeros);

      auto taylor3 = fun_threads[tid].Reverse(3, w);

    if(dense) {
    // store dense
      for(int Dj=0; Dj<Nparams; Dj++){
        int indexHere = 3*Dj;
//        double secondHere = taylor3[indexHere+1];
        gradientOut(Dj, Dk) = taylor3[indexHere+2];
        hessianOut(Dj, Dk) = taylor3[indexHere+1];
        thirdDiagOut(Dj, Dk) = taylor3[indexHere];//2*(taylor3[indexHere] - secondHere);
      }
    } else { // not dense

      // fill in the T_kii
      // first column is T_kii/2 + H_ki
      const int diagStart = DiagP[Dk];
      const int diagEnd = DiagP[Dk+1];

      for(int Di=diagStart; Di<diagEnd; Di++){
        int indexHere = 3*DiagRow[Di];
        double secondHere = taylor3[indexHere+1];
        hessianOut[Di] = secondHere;
        thirdDiagOut[Di] = taylor3[indexHere];//2*(taylor3[indexHere] - secondHere);
      } // Di
    } // else is sparse

  } // for k diagonal bit

#ifdef DEBUG
 if (verbose ) Rcpp::Rcout << "t" << tid;
#endif

} // parallel

if (verbose ) {
  Rcpp::Rcout << "done\n";
}


Rcpp::List resultList = Rcpp::List::create(
  Rcpp::Named("first") = gradientOut,
  Rcpp::Named("second") = hessianOut,
  Rcpp::Named("diag") = thirdDiagOut
  );

return resultList;

}



//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix thirdNonDiagonalsSparsity(
  Rcpp::NumericVector parameters, // beta, gamma, theta
  Rcpp::List data, 
  Rcpp::List config,
  Rcpp::IntegerMatrix pairs
  ) {

  int num_threads = 1;
  if (config.containsElementNamed("num_threads"))
    num_threads = Rcpp::as<int>(config["num_threads"]);

  const int Npairs = pairs.nrow();
  const int Nparams = parameters.size();
  Rcpp::LogicalMatrix result(Npairs, Nparams);

  // set up autodiff function
  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);  
  CppAD::ADFun<double> fun(ad_params, y);

  std::vector<double> x_val(Nparams);
  std::vector<double> y_val(1);
  std::vector<double> w{0.0, 0.0, 1.0};  

  std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
  for (size_t i = 0; i < num_threads; ++i) {
    fun_threads[i] = fun;
  }

  for (size_t i = 0; i < Nparams; ++i) {
    x_val[i] = parameters[i];
  }


    omp_set_num_threads(num_threads);
  #pragma omp parallel
  {
    const int tid=omp_get_thread_num();
    fun_threads[tid].Forward(0, x_val);

      #pragma omp for
  for(size_t Dpair=0;Dpair<Npairs;Dpair++) {

  std::vector<double> direction1(Nparams, 0.0), direction2(Nparams, 0.0);
  direction1[pairs(Dpair,0)]  = direction2[pairs(Dpair, 1)] = 1.0;     

//  fun_threads[tid].Forward(0, x_val);
  fun_threads[tid].Forward(1, direction1);
  fun_threads[tid].Forward(2, direction2);
  auto taylor3 = fun_threads[tid].Reverse(3, w);
      for(size_t Dk = 0;Dk < Nparams;Dk++) {
        if (!CppAD::NearEqual(taylor3[3*Dk], 0.0, 1e-12, 1e-12)) {
        result(Dpair, Dk) = 1L;
      } else {
        result(Dpair, Dk) = 0L;
      }
    } // Dk
  }
} // parallel
  return result;

}

/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix thirdOffDiagonals(
  Rcpp::NumericVector parameters, // beta, gamma, theta
  Rcpp::List data, 
  Rcpp::List config) {

  int num_threads = 1;
  if (config.containsElementNamed("num_threads"))
    num_threads = Rcpp::as<int>(config["num_threads"]);

  bool verbose = false;
  if (config.containsElementNamed("verbose")) {
   verbose = Rcpp::as<bool>(config["verbose"]);
 }
 bool dense = false;
 if (config.containsElementNamed("dense")) {
   dense = Rcpp::as<bool>(config["dense"]);
 }

 const int Nparams = parameters.size();

// hessian of random effects, lower triangle only, column format

// sparsity for third deriv
 Rcpp::List sparsity = config["sparsity"];
 Rcpp::List thirdList = sparsity["third"];
 Rcpp::List pairs = thirdList["pairs"];     

 Rcpp::IntegerVector sparsityIjI = pairs["i"];
 Rcpp::IntegerVector sparsityIjJ = pairs["j"];
 const int Npairs = sparsityIjJ.size();

 // only used for sprse
 Rcpp::IntegerVector sparsityIjP, sparsityIjPend;
 Rcpp::List sparsityThirdIjk;// = thirdList["ijk"];   
 Rcpp::IntegerVector sparsityIjkK;// = sparsityThirdIjk["k"];


// output
 Rcpp::NumericMatrix Tijk;

 if(dense){
    Tijk = Rcpp::NumericMatrix(Nparams, Npairs);
  } else {
    sparsityThirdIjk = thirdList["ijk"]; 
    sparsityIjP = pairs["p"];
    sparsityIjPend = pairs["pEnd"];
    sparsityIjkK = sparsityThirdIjk["k"];
    Tijk = Rcpp::NumericMatrix(sparsityIjkK.size(), 1L);
  }

if (verbose ) Rcpp::Rcout << "objects allocated\n";

// set up autodiff function
  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);  
  CppAD::ADFun<double> fun(ad_params, y);
  std::vector<double> x_val(Nparams), y_val(1),  w{0.0, 0.0, 1.0};  
  for (size_t i = 0; i < Nparams; ++i) {
    x_val[i] = parameters[i];
  }
  fun.Forward(0, x_val);

  // Replicate fun object for each thread
  std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
  for (int i = 0; i < num_threads; ++i) {
    fun_threads[i] = fun;
  }
  if (verbose ) {
    Rcpp::Rcout << "starting 3rd deriv, " << num_threads << " threads " <<
    Npairs << " jk pairs " << Nparams << " parameters " << "\n";
  }
  omp_set_num_threads(num_threads);

  #pragma omp parallel
  {
    const int tid=omp_get_thread_num();
    fun_threads[tid].Forward(0, x_val);

    std::vector<double> direction1(Nparams, 0.0), direction2(Nparams, 0.0);

// off diag T_ijk, pair is ij
  #pragma omp for
 for(int Dpair=0; Dpair < Npairs; ++Dpair) {

  std::fill(direction1.begin(), direction1.end(), 0.0);
  std::fill(direction2.begin(), direction2.end(), 0.0);

  const int Dj = sparsityIjJ[Dpair];
  const int Di = sparsityIjI[Dpair];
  direction1[Di]  = 1.0;
  direction2[Dj] = 1.0;     

  fun_threads[tid].Forward(1, direction1);
  fun_threads[tid].Forward(2, direction2);
  auto taylor3 = fun_threads[tid].Reverse(3, w);  

  // first column is third deriv combination
  //  T_iik/2 + T_jjk/2 + T_ijk 
  // columns of diag are the double deriv
  // rows of taylor3 are i

  if(dense) {
    for(int Dk=0; Dk<Nparams; Dk++){
      Tijk(Dk, Dpair) = taylor3[3*Dk];
    }
  } else { // sparse
    int DinIjkStart = sparsityIjP[Dpair];
    int DinIjkEnd = sparsityIjPend[Dpair];
    for(int DinIjk=DinIjkStart; 
      DinIjk<DinIjkEnd; DinIjk++
      ){
        int Dk = sparsityIjkK[DinIjk];
        Tijk[DinIjk] = taylor3[3*Dk];
    }
  } // end sparse
 } // Dk

} // parallel

if (verbose ) {
  Rcpp::Rcout << "done\n";
}


return Tijk;

}
