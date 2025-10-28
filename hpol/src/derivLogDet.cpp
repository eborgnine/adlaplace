#ifdef UNDEF

#include"hpol.hpp"
#include"loglikHelpers.hpp"
#include<omp.h>

//#define DEBUG



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


CppAD::vector<CppAD::AD<double>> logLikNoQ(
  CppAD::vector<CppAD::AD<double>> ad_params, 
  const Data &data, 
  const Config &config) {

  auto latent=unpack_params(ad_params, data, config);

  Rcpp::NumericVector Sstrata(2);
  Sstrata[0]=0;
  Sstrata[1]=data.Nstrata;

  auto result = loglikSeq(Sstrata, latent.gamma, latent, data, config);
  return(result);
};  



CppAD::vector<CppAD::AD<double>>  logLikQ(
    CppAD::vector<CppAD::AD<double>> ad_params, 
    const Data& data,
    const Config& config
) {

  auto latent = unpack_params(ad_params, data, config);

    CppAD::vector<CppAD::AD<double>> gammaScaled(data.Ngamma);
    CppAD::vector<CppAD::AD<double>> result(1);
    result[0] = 0;

    for (size_t D = 0; D < data.Ngamma; ++D) {
        size_t mapHere = data.map[D];
        CppAD::AD<double> thetaHere = latent.theta[mapHere];
        CppAD::AD<double> logThetaHere = latent.logTheta[mapHere];

        gammaScaled[D] = latent.gamma[D] / thetaHere;

        result[0] += logThetaHere +
                      (0.5 * data.Qdiag[D]) * gammaScaled[D] * gammaScaled[D] ;
    }

      // Q offdiag    
    for(size_t D = 0; D < data.Nq; D++) {
        result[0] += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] 
          * (data.QsansDiag.x[D]);
    }


  return(result);
}


template<class Type>
CppAD::vector<Type>  objectiveFunctionInternal(
 const CppAD::vector<Type>& ad_params,  
 const Data& data,
 const Config& config
 );



/* function for third derivatives, of laplace approx */

Rcpp::List thirdDiagonals(
  Rcpp::NumericVector parameters, // beta, gamma, theta
  Rcpp::List dataList, 
  Rcpp::List configList // sparsity$second full and random sparsity$third jk and ijk
  ) {

  Data   data(dataList);
  Config config(configList);

  const int Nparams = parameters.size();

// output
// Rcpp::NumericVector gradient(Nparams);
  Rcpp::NumericMatrix hessianOut, thirdDiagOut;
  Rcpp::List hessianList, hessian, diagSparsity;
  Rcpp::IntegerVector Hrow, Hp,  DiagRow, DiagP;
  Rcpp::List sparsity = config.sparsity; 

  const std::vector<double> x_val(parameters.begin(), parameters.end());

  bool dense = config.dense;
  if(!sparsity.size()) {
    if(!dense) {
      Rcpp::warning("no sparsity provided, switching to dense\n");
      dense = true;
    }
  }

  Rcpp::NumericVector gradientOut(Nparams);

  if(dense){
    hessianOut = Rcpp::NumericMatrix(Nparams, Nparams);
    thirdDiagOut = Rcpp::NumericMatrix(Nparams, Nparams);
  } else {
// hessian of random effects, lower triangle only, column format
    hessianList = sparsity["second"];
    hessian = hessianList["full"];
    diagSparsity = hessianList["nonSymmetric"];
    Hrow = hessian["i"]; 
    Hp = hessian["p"];
    DiagRow = diagSparsity["i"]; 
    DiagP = diagSparsity["p"];

    hessianOut = Rcpp::NumericMatrix(DiagRow.size(), 1);
    thirdDiagOut = Rcpp::NumericMatrix(DiagRow.size(), 1);
  }

  if (config.verbose ) Rcpp::Rcout << "objects allocated\n";

//


  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );


  if (config.verbose ) Rcpp::Rcout << "starting parallel " << config.num_threads << " threads\n";

  #pragma omp parallel
  {
 // set up autodiff function

//        const int tid=omp_get_thread_num();


    std::vector<double> y_val(1);
    std::vector<double> w{0.0, 0.0, 1.0};  
    std::vector<double> direction(Nparams, 0.0), directionZeros(Nparams, 0.0);


    CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
    for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = x_val[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal<CppAD::AD<double>>(ad_params, data, config);  

  CppAD::ADFun<double> fun(ad_params, y);

  y_val = fun.Forward(0, x_val);


  // diagonal entries to subtract off
  // computing T_kii and H_ki, columnns are i, rows are k
    #pragma omp for
  for(int Dk=0; Dk < Nparams; ++Dk) {

    std::fill(direction.begin(), direction.end(), 0.0);
    direction[Dk]  = 1.0;     

    fun.Forward(1, direction);
    fun.Forward(2, directionZeros);

    auto taylor3 = fun.Reverse(3, w);

    if(dense) {
    // store dense
      for(int Dj=0; Dj<Nparams; 
        Dj++){
        int indexHere = 3*Dj;
      hessianOut(Dj, Dk) = taylor3[indexHere+1];
      thirdDiagOut(Dj, Dk) = taylor3[indexHere];
    }

    } else { // not dense

      // fill in the T_kii
      // first column is T_kii/2 + H_ki
      const int diagStart = DiagP[Dk];
      const int diagEnd = DiagP[Dk+1];

      for(int Di=diagStart; Di<diagEnd; Di++){
        int indexHere = 3*DiagRow[Di];
        hessianOut[Di] = taylor3[indexHere+1];
        thirdDiagOut[Di] = taylor3[indexHere];
      } // Di
    } // else is sparse

    if(Dk == 0) {
      for(int Dj=0; Dj<Nparams; Dj++){
        int indexHere = 3*Dj;
        gradientOut[Dj] = taylor3[indexHere+2];
      }
    }
  } // for k diagonal bit

} // parallel

if (config.verbose ) {
  Rcpp::Rcout << "done\n";
}


Rcpp::List resultList = Rcpp::List::create(
  Rcpp::Named("first") = gradientOut,
  Rcpp::Named("second") = hessianOut,
  Rcpp::Named("diag") = thirdDiagOut
  );

return resultList;

}




Rcpp::LogicalMatrix thirdNonDiagonalsSparsity(
  Rcpp::NumericVector parameters, // beta, gamma, theta
  Rcpp::List dataList, 
  Rcpp::List configList,
  Rcpp::IntegerMatrix pairs
  ) {

  Data   data(dataList);
  Config config(configList);


  const int Npairs = pairs.nrow();
  const int Nparams = parameters.size();

const Rcpp::IntegerVector colI = pairs.column(0);
const Rcpp::IntegerVector colJ = pairs.column(1);

const std::vector<int> pairsI(colI.begin(), colI.end());
const std::vector<int> pairsJ(colJ.begin(), colJ.end());

  Rcpp::LogicalMatrix result(Npairs, Nparams);

  omp_set_num_threads(config.num_threads);
  #pragma omp parallel
  {

  // set up autodiff function
  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
  std::vector<double> x_val(Nparams);
  std::vector<double> y_val(1);

  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
    x_val[D] = parameters[D];
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  auto y = logLikNoQ(ad_params, data, config);  
  CppAD::ADFun<double> fun(ad_params, y);
  fun.Forward(0, x_val);


  const std::vector<double> direction2(Nparams, 0.0);
  const std::vector<double> w{0.0, 0.0, 1.0};  

      #pragma omp for
    for(size_t Dpair=0;Dpair<Npairs;Dpair++) {

      std::vector<double> direction1(Nparams, 0.0);
      direction1[pairsI[Dpair] ]  = direction1[pairsJ[Dpair] ] = 1.0;     

//  fun_threads[tid].Forward(0, x_val);
      fun.Forward(1, direction1);
      fun.Forward(2, direction2);
      auto taylor3 = fun.Reverse(3, w);
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

Rcpp::NumericMatrix thirdOffDiagonals(
  Rcpp::NumericVector parameters, // beta, gamma, theta
  Rcpp::List dataList, 
  Rcpp::List configList) {



  Data   data(dataList);
  Config config(configList);


  const int Nparams = parameters.size();

// hessian of random effects, lower triangle only, column format

// sparsity for third deriv
  const Rcpp::List sparsity = config.sparsity;
  const Rcpp::List thirdList = sparsity["third"];
  const Rcpp::List pairs = thirdList["pairs"];     

  const Rcpp::IntegerVector sparsityIjI = pairs["i"];
  const Rcpp::IntegerVector sparsityIjJ = pairs["j"];
  const Rcpp::LogicalVector pairsNoData = pairs["nodata"];
  const Rcpp::LogicalVector pairsNoQ = pairs["noQ"];
  const int Npairs = sparsityIjJ.size();

 // only used for sparse
  Rcpp::IntegerVector sparsityIjP, sparsityIjPend;
  Rcpp::List sparsityThirdIjk;// = thirdList["ijk"];   
  Rcpp::IntegerVector sparsityIjkK;// = sparsityThirdIjk["k"];

// output
  Rcpp::NumericMatrix Tijk;

  if(config.dense){
    Tijk = Rcpp::NumericMatrix(Nparams, Npairs);
    if (config.verbose ) Rcpp::Rcout << "Nparams " << Nparams << " Npairs " << Npairs << "\n";
  } else {
    sparsityThirdIjk = thirdList["ijk"]; 
    sparsityIjP = pairs["p"];
    sparsityIjPend = pairs["pEnd"];
    sparsityIjkK = sparsityThirdIjk["k"];
    Tijk = Rcpp::NumericMatrix(sparsityIjkK.size(), 2L);
  }

  if (config.verbose ) {
    Rcpp::Rcout << "starting 3rd deriv, " << config.num_threads << " threads " <<
    Npairs << " jk pairs " << Nparams << " parameters " << "\n";
  }

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
// set up autodiff function

  CppAD::vector<CppAD::AD<double>> ad_params(Nparams), ad_params_Q(Nparams);
  std::vector<double> x_val(Nparams);  

  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
    ad_params_Q[D] = parameters[D];  // Initialize CppAD variables
    x_val[D] = parameters[D];
  }

  CppAD::Independent(ad_params);  
  auto y = logLikNoQ(ad_params, data, config);  
  CppAD::ADFun<double> fun(ad_params, y);
  fun.Forward(0, x_val);

  CppAD::Independent(ad_params_Q);  
  auto yQ = logLikQ(ad_params_Q, data, config);
  CppAD::ADFun<double> funQ(ad_params_Q, yQ);
  funQ.Forward(0, x_val);


    const std::vector<double>  w{0.0, 0.0, 1.0};  
    const std::vector<double> direction2(Nparams, 0.0);
    std::vector<double> direction1(Nparams, 0.0);


// off diag T_ijk, pair is ij
  #pragma omp for
    for(int Dpair=0; Dpair < Npairs; ++Dpair) {

    std::vector<double> taylor3, taylor3Q;

        const int Di = sparsityIjI[Dpair];
        const int Dj = sparsityIjJ[Dpair];
        const bool havedata = !pairsNoData[Dpair];
        const bool haveQ = !pairsNoQ[Dpair];

        std::fill(direction1.begin(), direction1.end(), 0.0);
        direction1[Di] = direction1[Dj] = 1.0;     

        if(havedata) {
          fun.Forward(1, direction1);
          fun.Forward(2, direction2);
          taylor3 = fun.Reverse(3, w);  
        }

        if(haveQ) {
          funQ.Forward(1, direction1);
          funQ.Forward(2, direction2);
          taylor3Q = funQ.Reverse(3, w);  
        }

  // first column is third deriv combination
  //  T_iik/2 + T_jjk/2 + T_ijk 
  // columns of diag are the double deriv
  // rows of taylor3 are i
    if(config.dense) {
      for(int Dk=0; Dk<Nparams; Dk++){
        if(havedata) {
          Tijk(Dk, Dpair) += taylor3[3*Dk];
        } else {
          Tijk(Dk, Dpair) =0;
        }
        if(haveQ) Tijk(Dk, Dpair) += taylor3Q[3*Dk];
      }
    } else { // sparse
      int DinIjkStart = sparsityIjP[Dpair];
      int DinIjkEnd = sparsityIjPend[Dpair];
      for(int DinIjk=DinIjkStart; 
          DinIjk<DinIjkEnd; DinIjk++
      ){
        int Dk = sparsityIjkK[DinIjk];
        if(havedata) {
          Tijk(DinIjk,0) = taylor3[3*Dk];
        } else {
          Tijk(DinIjk,0) = 0;
        }
        if(haveQ) {
          Tijk(DinIjk,1) = taylor3Q[3*Dk];
        } else {
          Tijk(DinIjk,1) = 0;
        }
      }
    } // end sparse
  } // Dpair

} // parallel


if (config.verbose ) {
  Rcpp::Rcout << "done\n";
}

return Tijk;

}
#endif
