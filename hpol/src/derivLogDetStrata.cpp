#include"hpol.hpp"
#include"loglikHelpers.hpp"
#include<omp.h>

//#define DEBUG





/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::List thirdDiagonalsStrata(
  const Rcpp::NumericVector parameters, // beta, gamma, theta
  const Rcpp::List data, 
  const Rcpp::List config, 
  const Rcpp::List sparsity,
  const Rcpp::List strata
  ) {

  const Data   dataC(data);
  const Config configC(config);

  const int Nparams = parameters.size();
  const int Ngroup = strata.size();

  const Rcpp::NumericVector strataI = strata["i"], strataP = strata["p"];


  bool dense = configC.dense;
  if(!sparsity.size()) {
    if(!dense) {
      Rcpp::warning("no sparsity provided, switching to dense\n");
      dense = true;
    }
  }
  Rcpp::List hessianOut(Ngroup), thirdDiagOut(Ngroup);

  for(size_t D=0;D<Ngroup;D++) {
    if(dense) {
      hessianOut[D] = Rcpp::NumericMatrix(Nparams, Nparams);
      thirdDiagOut[D] = Rcpp::NumericMatrix(Nparams, Nparams);
    } else {
      Rcpp::List sparseHere = sparsity[D];
      Rcpp::List twoHere = sparseHere["two"];
      Rcpp::List fullHere = twoHere["full"];
      Rcpp::IntegerVector fullIhere = fullHere["i"];

      const size_t NoutRowsH = fullIhere.size();

      Rcpp::List twoThere = twoHere["nonSymmetric"];
      Rcpp::List twoTihere = twoThere["i"];

      const size_t NoutRowsT = twoTihere.size();

      hessianOut[D] = Rcpp::NumericMatrix(NoutRowsH, 1L);
      thirdDiagOut[D] = Rcpp::NumericMatrix(NoutRowsT, 1L);
    }

  }
  Rcpp::NumericMatrix gradientOut(Nparams, Ngroup);

  const std::vector<double> x_val(parameters.begin(), parameters.end());


  if (configC.verbose ) Rcpp::Rcout << "objects allocated\n";

  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );


  if (configC.verbose ) Rcpp::Rcout << "starting parallel " << configC.num_threads << " threads\n";

  #pragma omp parallel
  {
 // set up autodiff function

    std::vector<double> y_val(1);
    const std::vector<double> w{0.0, 0.0, 1.0};  
    std::vector<double> direction(Nparams, 0.0), directionZeros(Nparams, 0.0);


    CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
    for (size_t D = 0; D < Nparams; D++) {
      ad_params[D] = x_val[D];  // Initialize CppAD variables
    }

  #pragma omp for
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      const Rcpp::List sparsityHere = sparsity[Dgroup];
      const Rcpp::List sparsitySecond = sparsityHere["second"];

      // DiagP = second$nonSymmetric$p  
      const Rcpp::List nonSymmetricHere = sparsitySecond["nonSymmetric"];
      const Rcpp::IntegerVector diagP = nonSymmetricHere["p"];
      const Rcpp::IntegerVector diagRow = nonSymmetricHere["i"];
      const Rcpp::List fullHere = sparsitySecond["full"];
      const Rcpp::IntegerVector hessianP = fullHere["p"];
      const Rcpp::IntegerVector hessianRow = fullHere["i"];


      Rcpp::NumericMatrix hessianOutHere = hessianOut[Dgroup];
      Rcpp::NumericMatrix thirdDiagOutHere = thirdDiagOut[Dgroup];


      CppAD::Independent(ad_params);  

      auto latent=unpack_params(ad_params, dataC, configC);

      CppAD::AD<double> loglik = 0;
      CppAD::vector<CppAD::AD<double>> minusLogDens(1);

      const size_t end = strataP[Dgroup+1];
      for (size_t Dindex = strataP[Dgroup]; Dindex < end;  Dindex++) {

        const size_t Dstrata = strata[Dindex];

        auto etaHere = compute_eta_for_stratum(
          Dstrata, dataC, latent.gamma, latent.beta);

        auto contrib = accumulate_contrib_for_stratum(
          Dstrata, dataC, etaHere, latent, configC);

        loglik += contrib[0];
    } // Dstrata

    minusLogDens[0] = -loglik;

    CppAD::ADFun<double> fun(ad_params, minusLogDens);

    y_val = fun.Forward(0, x_val);

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
        hessianOutHere(Dj, Dk) = taylor3[indexHere+1];
        thirdDiagOutHere(Dj, Dk) = taylor3[indexHere];
      }
    } else { // not dense

      // fill in the T_kii
      // first column is T_kii/2 + H_ki

      const int diagEnd = diagP[Dk+1];

      for(int Di=diagP[Dk]; Di<diagEnd; Di++){
        int indexHere = 3*diagRow[Di];
        thirdDiagOutHere[Di] = taylor3[indexHere+1];
      }
      const int hessianEnd = hessianP[Dk+1];
      for(int Di=hessianP[Dk]; Di<hessianEnd; Di++){
        int indexHere = 3*hessianRow[Di];
        hessianOutHere[Di] = taylor3[indexHere];
      } // Di
    } // else is sparse

    if(Dk == 0) {
      for(int Dj=0; Dj<Nparams; Dj++){
        int indexHere = 3*Dj;
        gradientOut(Dj,Dgroup) = taylor3[indexHere+2];
      }
    }
  } // for k diagonal bit
  } //group

} // parallel

if (configC.verbose ) {
  Rcpp::Rcout << "done\n";
}


Rcpp::List resultList = Rcpp::List::create(
  Rcpp::Named("first") = gradientOut,
  Rcpp::Named("second") = hessianOut,
  Rcpp::Named("diag") = thirdDiagOut
  );

return resultList;

}




/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::List thirdOffDiagonalsStrata(
  const Rcpp::NumericVector parameters, // beta, gamma, theta
  const Rcpp::List data, 
  const Rcpp::List config,
  const Rcpp::List sparsity, 
  const Rcpp::List strata) {

  const Data   dataC(data);
  const Config configC(config);


  const int Nparams = parameters.size();
  const int Ngroups = strata.size();
  const Rcpp::NumericVector strataI = strata["i"], strataP = strata["p"];

// output
  Rcpp::List Tijk(Ngroups);

  for(size_t D=0;D<Ngroups;D++) {

    const Rcpp::List sparsityHere = sparsity[D];
    const Rcpp::List threeHere = sparsityHere["three"];
    const Rcpp::List pairsHere = threeHere["pairs"];
    const Rcpp::IntegerVector pairsI = pairsHere["i"];
    const size_t Npairs = pairsI.size();

    if(configC.dense) {
      Tijk[D] = Rcpp::NumericMatrix(Nparams, Npairs);
    } else {
      const Rcpp::List ijkHere = threeHere["ijk"];
      const Rcpp::List sparsityIjkK = ijkHere["k"];
      Tijk[D] = Rcpp::NumericMatrix(sparsityIjkK.size(), 1L); 
    }
  }

  const std::vector<double> x_val(parameters.begin(), parameters.end());


  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
// set up autodiff function

  const std::vector<double>  w{0.0, 0.0, 1.0};  
  const std::vector<double> direction2(Nparams, 0.0);
  std::vector<double> direction1(Nparams, 0.0);

    CppAD::vector<CppAD::AD<double>> ad_params(Nparams);

    for (size_t D = 0; D < Nparams; D++) {
      ad_params[D] = x_val[D];
    }

  #pragma omp for
    for(size_t Dgroup = 0; Dgroup < Ngroups; ++Dgroup) {

      const Rcpp::List sparsityHere = sparsity[Dgroup];
      const Rcpp::List threeHere = sparsityHere["three"];
      const Rcpp::List pairsHere = threeHere["pairs"];
      const Rcpp::List ijkHere = threeHere["ijk"];
      const Rcpp::IntegerVector pairsI = pairsHere["i"];
      const Rcpp::IntegerVector pairsJ = pairsHere["j"];
      const Rcpp::IntegerVector pairsP = pairsHere["p"];
      const Rcpp::IntegerVector pairsPend = pairsHere["pEnd"];
      const Rcpp::IntegerVector sparsityIjk = ijkHere["k"];


      const size_t Npairs = pairsI.size();

      Rcpp::NumericMatrix TijkHere = Tijk[Dgroup];

      CppAD::AD<double> loglik = 0;
      CppAD::vector<CppAD::AD<double>> minusLogDens(1);

      CppAD::Independent(ad_params);
      auto latent=unpack_params(ad_params, dataC, configC);

      const size_t end = strataP[Dgroup+1];
      for (size_t Dindex = strataP[Dgroup]; Dindex < end;  Dindex++) {

        const size_t Dstrata = strata[Dindex];

        auto etaHere = compute_eta_for_stratum(
          Dstrata, dataC, latent.gamma, latent.beta);

        auto contrib = accumulate_contrib_for_stratum(
          Dstrata, dataC, etaHere, latent, configC);

        loglik += contrib[0];
    } // Dstrata

    minusLogDens[0] = -loglik;

  CppAD::ADFun<double> fun(ad_params, minusLogDens);
  fun.Forward(0, x_val);

  for(int Dpair=0; Dpair < Npairs; ++Dpair) {

    const int Di = pairsI[Dpair];
    const int Dj = pairsJ[Dpair];

    std::fill(direction1.begin(), direction1.end(), 0.0);
    direction1[Di] = direction1[Dj] = 1.0;     

      fun.Forward(1, direction1);
      fun.Forward(2, direction2);
      auto taylor3 = fun.Reverse(3, w);  


  // first column is third deriv combination
  //  T_iik/2 + T_jjk/2 + T_ijk 
  // columns of diag are the double deriv
  // rows of taylor3 are i
    if(configC.dense) {
      for(int Dk=0; Dk<Nparams; Dk++){
          TijkHere(Dk, Dpair) += taylor3[3*Dk];
      }
    } else { // sparse
      int DinIjkStart = pairsP[Dpair];
      int DinIjkEnd = pairsPend[Dpair];
      for(int DinIjk=DinIjkStart; DinIjk<DinIjkEnd; DinIjk++){
        int Dk = sparsityIjk[DinIjk];
        TijkHere(DinIjk,0) = taylor3[3*Dk];
      }
    } // end sparse
  } // Dpair
} //Dgroup
} // parallel



return Tijk;

}
