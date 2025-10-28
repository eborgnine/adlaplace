#include"hpol.hpp"



/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::List thirdDiagonalsStrata(
  const Rcpp::NumericVector parameters, // beta, gamma, theta
  const Rcpp::List data, 
  const Rcpp::List config
  ) {

  const Data   dataC(data);
  const Config configC(config);

  const Rcpp::List sparsity = configC.group_sparsity;
  const Rcpp::List strata = configC.groups;

  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];
  const int Ngroup = strataP.size()-1;

struct ThirdPack {
  std::vector<int> pairsI, pairsJ, pairsP, pairsPend, ijkK;
};
std::vector<ThirdPack> third(Ngroup);

for (int g = 0; g < Ngroup; ++g) {
  Rcpp::List spHere   = sparsity[g];
  Rcpp::List three    = spHere["third"];
  Rcpp::List pairs    = three["pairs"];
  Rcpp::List ijk      = three["ijk"];

  third[g].pairsI   = Rcpp::as<std::vector<int>>(pairs["i"]);
  third[g].pairsJ   = Rcpp::as<std::vector<int>>(pairs["j"]);
  third[g].pairsP   = Rcpp::as<std::vector<int>>(pairs["p"]);
  third[g].pairsPend= Rcpp::as<std::vector<int>>(pairs["pEnd"]);
  third[g].ijkK     = Rcpp::as<std::vector<int>>(ijk["k"]);
}

  const int Nparams = parameters.size();
  const int NparamsSq = Nparams*Nparams;
  if(Nparams == dataC.Ngamma) Rcpp::warning("need full parameters, only random has been sent");

  bool dense = configC.dense;
  const std::vector<double> x_val(parameters.begin(), parameters.end());

  auto fun = getAdFun(x_val, dataC, configC);
  auto Qfun = getAdFunQ(x_val, dataC, configC); 

  const Rcpp::List secondAll = configC.sparsity["second"];
  const Rcpp::List fullAll = secondAll["full"];
  const Rcpp::List nsAll = secondAll["nonSymmetric"];

  const Rcpp::IntegerVector iFullAll = fullAll["i"];
  const Rcpp::IntegerVector iNsAll = nsAll["i"];


  const size_t NoutRowsT = dense?NparamsSq:iNsAll.size();
  const size_t NoutRowsH = dense?NparamsSq:iFullAll.size();

  std::vector<double> hessianOut(NoutRowsH, 0.0);
  std::vector<double> thirdDiagOut(NoutRowsT, 0.0);
  std::vector<double> gradientOut(Nparams, 0.0);
  std::vector<std::vector<double>> TijkOut(Ngroup);

  if (configC.verbose ) Rcpp::Rcout << "starting parallel " << configC.num_threads << " threads\n";

  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {

    std::vector<double> hessianOutHere(NoutRowsH, 0.0);
    std::vector<double> thirdDiagOutHere(NoutRowsT, 0.0);
    std::vector<double> gradientOutHere(Nparams, 0.0);
    std::vector<double> denseThirdHere(NparamsSq);

    const std::vector<double> w{0.0, 0.0, 1.0};  
    std::vector<double> direction(Nparams, 0.0);
    const std::vector<double> directionZeros(Nparams, 0.0);

# pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      fun[Dgroup].fun.Forward(0, x_val);

      auto& iNS = fun[Dgroup].nsRowCol[0];
      auto& matchNS = fun[Dgroup].nsRowCol[2];
      auto& pNS = fun[Dgroup].nsP;

      auto& iFull = fun[Dgroup].outRowCol[0];
      auto& matchFull = fun[Dgroup].outRowCol[2];
      auto& pFull = fun[Dgroup].outP;

//      Rcpp::Rcout << "group " << Dgroup << " " << iNS[0] << "a " << matchNS[0] << "b " << pNS[0] << "c " << iFull[0] << "d " << matchFull[0] << "e " << pFull[0] << "kk" << pFull[1] << "f\n";

      for(int Dk=0; Dk < Nparams; ++Dk) {

        std::fill(direction.begin(), direction.end(), 0.0);
        direction[Dk]  = 1.0;     

        fun[Dgroup].fun.Forward(1, direction);
        fun[Dgroup].fun.Forward(2, directionZeros);

        auto taylor3 = fun[Dgroup].fun.Reverse(3, w);

    // store dense
        const size_t colStart = Dk * Nparams;
        for(int Dj=0; Dj<Nparams; Dj++){
          const size_t indexHere = 3*Dj;
          const size_t rowHere = colStart + Dj;
          // first column is T_kii/2 + H_ki
          denseThirdHere[rowHere] = 2*(taylor3[indexHere] - taylor3[indexHere+1]);
          if(dense) {
            hessianOutHere[rowHere] += taylor3[indexHere+1];
            thirdDiagOutHere[rowHere] += denseThirdHere[rowHere];
          }
        }
      if(!dense) { // not dense
      // fill in the T_kii
      // first column is T_kii/2 + H_ki
        const size_t nsEnd = pNS[Dk+1];
        for(size_t Di=pNS[Dk]; Di<nsEnd; Di++){
          const size_t indexHere = 3*iNS[Di];
          const size_t allHere = matchNS[Di];
          thirdDiagOutHere[allHere] += 2*(taylor3[indexHere] - taylor3[indexHere+1]);
        }
        const size_t fullEnd = pFull[Dk+1];
 
        for(size_t Di=pFull[Dk]; Di<fullEnd; Di++){
          const size_t indexHere = 3*iFull[Di];
          const size_t allHere = matchFull[Di];
          hessianOutHere[allHere] += taylor3[indexHere+1];
        } // Di
      } // else is sparse

      if(Dk == 0) {
        for(int Dj=0; Dj<Nparams; Dj++){
          const int indexHere = 3*Dj;
          gradientOutHere[Dj] += taylor3[indexHere+2];
        }
      }
    } // for k diagonal bit

    // off diagonals

      auto& pairsI = third[Dgroup].pairsI;
      auto& pairsJ = third[Dgroup].pairsJ;
      auto& pairsP = third[Dgroup].pairsP;
      auto& pairsPend = third[Dgroup].pairsPend;
      auto& sparsityIjk = third[Dgroup].ijkK;

      const size_t Npairs = pairsI.size();

      const size_t NtijkHere = dense?Nparams*Npairs:sparsityIjk.size();
      std::vector<double> TijkHere(NtijkHere);

      for(int Dpair=0; Dpair < Npairs; ++Dpair) {

        const size_t Di = pairsI[Dpair];
        const size_t Dj = pairsJ[Dpair];
        const size_t pairStart = pairsP[Dpair];
        const size_t DiNparams = Di* Nparams;
        const size_t DjNparams = Dj* Nparams;
        const size_t DcolHere =  Dpair * Nparams;

        std::fill(direction.begin(), direction.end(), 0.0);
        direction[Di] = direction[Dj] = 1.0;     

        fun[Dgroup].fun.Forward(1, direction);
        fun[Dgroup].fun.Forward(2, directionZeros);
        auto taylor3 = fun[Dgroup].fun.Reverse(3, w);

  // first column is third deriv combination
      //  T_iik/2 + T_jjk/2 + T_ijk 

        if(dense) {
          for(size_t Dk=0; Dk<Nparams; Dk++){
            const double TiikTjjk = thirdDiagOutHere[DiNparams + Dk] + thirdDiagOutHere[DjNparams + Dk];
            TijkHere[DcolHere + Dk] = taylor3[3*Dk] - TiikTjjk/2; 
          }
        } else {
          const size_t NthisPair = pairsPend[Dpair] - pairsP[Dpair];
          for(size_t Dindex=0; Dindex<NthisPair; Dindex++){
            const size_t DindexInIjk = pairStart + Dindex;
            const size_t Dk = sparsityIjk[DindexInIjk];
            const double TiikTjjk = thirdDiagOutHere[DiNparams + Dk] + thirdDiagOutHere[DjNparams + Dk];
            TijkHere[pairStart + Dindex] = taylor3[3*Dk] - TiikTjjk/2;
          }           
        }
      } // Dpair
        // save the TijkHere
      TijkOut[Dgroup] = TijkHere;
  } //group


# pragma omp single
  {
        Qfun.fun.Forward(0, x_val);

        auto& iNS = Qfun.nsRowCol[0];
        auto& matchNS = Qfun.nsRowCol[2];
        auto& pNS = Qfun.nsP;
        auto& iFull = Qfun.outRowCol[0];
        auto& matchFull = Qfun.outRowCol[2];
        auto& pFull = Qfun.outP;


      for(int Dk=0; Dk < Nparams; ++Dk) {
        std::fill(direction.begin(), direction.end(), 0.0);
        direction[Dk]  = 1.0;     
        Qfun.fun.Forward(1, direction);
        Qfun.fun.Forward(2, directionZeros);

        auto taylor3 = Qfun.fun.Reverse(3, w);


      if(dense) {
    // store dense
        const size_t colStart = Dk * Nparams;
        for(int Dj=0; Dj<Nparams; Dj++){
          const size_t indexHere = 3*Dj;
          const size_t rowHere = colStart + Dj;

          hessianOutHere[rowHere] += taylor3[indexHere+1];
          thirdDiagOutHere[rowHere] += taylor3[indexHere];
        }
      } else { // not dense

        const size_t nsEnd = pNS[Dk+1];
        for(size_t Di=pNS[Dk]; Di<nsEnd; Di++){
          const size_t indexHere = 3*iNS[Di];
          const size_t allHere = matchNS[Di];
          thirdDiagOutHere[allHere] += taylor3[indexHere];
        } // Di

        const size_t fullEnd = pFull[Dk+1];
        for(size_t Di=pFull[Dk]; Di<fullEnd; Di++){
          const size_t indexHere = 3*iFull[Di];
          const size_t allHere = matchFull[Di];
          hessianOutHere[allHere] += taylor3[indexHere+1];
        }
      } // else is sparse
    } // Dk
  } // Q

# pragma omp critical
  {
    for(size_t Dcol=0;Dcol<Nparams;Dcol++) {
      gradientOut[Dcol] += gradientOutHere[Dcol];
    }
    for(size_t D=0;D<NoutRowsH;D++) {
      hessianOut[D] += hessianOutHere[D];
    }
    for(size_t D=0;D<NoutRowsT;D++) {
      thirdDiagOut[D] += thirdDiagOutHere[D];
    }
  } // critical
} // parallel


Rcpp::List TijkR(Ngroup);
for(size_t D=0;D<Ngroup;D++) {
  TijkR[D] = Rcpp::NumericVector(TijkOut[D].begin(), TijkOut[D].end());
}

if (configC.verbose ) {
  Rcpp::Rcout << "done\n";
}

Rcpp::NumericVector gradientOutR(gradientOut.begin(), gradientOut.end());
Rcpp::NumericVector hessianOutR(hessianOut.begin(), hessianOut.end());
Rcpp::NumericVector thirdDiagOutR(thirdDiagOut.begin(), thirdDiagOut.end());

Rcpp::List resultList = Rcpp::List::create(
  Rcpp::Named("first") = gradientOutR,
  Rcpp::Named("second") = hessianOutR,
  Rcpp::Named("diag") = thirdDiagOutR,
  Rcpp::Named("Tijk") = TijkR
  );


return resultList;

}


#ifdef UNDEF

/* function for third derivatives, of laplace approx 
 ' @export
 [[Rcpp::export]]*/
Rcpp::List thirdOffDiagonalsStrata(
  const Rcpp::NumericVector parameters, // beta, gamma, theta
  const Rcpp::List data, 
  const Rcpp::List config) {

  const Data   dataC(data);
  const Config configC(config);

  const Rcpp::List sparsity = configC.group_sparsity;
  const Rcpp::List strata = configC.groups;


  const int Nparams = parameters.size();
  const int Ngroups = strata.size();
  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

// output
  Rcpp::List Tijk(Ngroups);

  for(size_t D=0;D<Ngroups;D++) {

    const Rcpp::List sparsityHere = sparsity[D];
    const Rcpp::List threeHere = sparsityHere["third"];
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

  #pragma omp for
    for(size_t Dgroup = 0; Dgroup < Ngroups; ++Dgroup) {

      CppAD::ADFun<double> fun = adFunGroup(x_val,  
        dataC, configC, strataI,
        strataP[Dgroup], strataP[Dgroup+1]);

      const Rcpp::List sparsityHere = sparsity[Dgroup];
      const Rcpp::List threeHere = sparsityHere["third"];
      const Rcpp::List pairsHere = threeHere["pairs"];
      const Rcpp::List ijkHere = threeHere["ijk"];
      const Rcpp::IntegerVector pairsI = pairsHere["i"];
      const Rcpp::IntegerVector pairsJ = pairsHere["j"];
      const Rcpp::IntegerVector pairsP = pairsHere["p"];
      const Rcpp::IntegerVector pairsPend = pairsHere["pEnd"];
      const Rcpp::IntegerVector sparsityIjk = ijkHere["k"];


      const size_t Npairs = pairsI.size();

      Rcpp::NumericMatrix TijkHere = Tijk[Dgroup];



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

// combine over groups

return Tijk;

}

#endif
