#include"hpol.hpp"


struct ThirdPack {
  std::vector<int> pairsI, pairsJ, pairsP, pairsPend, ijkK;
};


Rcpp::List thirdStrata(
  const std::vector<double>&  parameters, // beta, gamma, theta
  const std::vector<ThirdPack>& third, 
  const Config& config,
  std::vector<GroupPack>& fun,
  GroupPack& Qfun
  ) {

  const std::vector<double>& x_val = parameters;

  const size_t Nparams = parameters.size();
  const size_t NparamsSq = Nparams*Nparams;
  const size_t Ngroup = fun.size();

  bool dense = config.dense;


  const Rcpp::List secondAll = config.sparsity["second"];
  const Rcpp::List fullAll = secondAll["full"];
  const Rcpp::List nsAll = secondAll["nonSymmetric"];
  const Rcpp::IntegerVector iFullAll = fullAll["i"];
  const Rcpp::IntegerVector iNsAll = nsAll["i"];
  const size_t NoutRowsT = dense?NparamsSq:iNsAll.size();
  const size_t NoutRowsH = dense?NparamsSq:iFullAll.size();
  const Rcpp::List thirdAll = config.sparsity["third"];
  const Rcpp::List pairsAll = thirdAll["pairs"];
  const Rcpp::NumericVector pairsIall = pairsAll["i"];
  const size_t NpairsAll = pairsIall.size();
  const size_t NpairsAllNparams = NpairsAll*Nparams;

  std::vector<double> hessianOut(NoutRowsH, 0.0);
  std::vector<double> thirdDiagOut(NoutRowsT, 0.0);
  std::vector<double> gradientOut(Nparams, 0.0);
  std::vector<std::vector<double>> TijkOut(Ngroup);

  if (config.verbose ) Rcpp::Rcout << "starting parallel " << config.num_threads << " threads\n";

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {

    std::vector<double> hessianOutHere(NoutRowsH, 0.0);
    std::vector<double> thirdDiagOutHere(NoutRowsT, 0.0);
    std::vector<double> gradientOutHere(Nparams, 0.0);

    const std::vector<double> w{0.0, 0.0, 1.0};  
    std::vector<double> direction(Nparams, 0.0);
    const std::vector<double> directionZeros(Nparams, 0.0);

    std::vector<double> TijkHere;


# pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      // keep all the third Tiik for this group
      std::vector<double> denseThirdHere(NparamsSq, 0.0);

      fun[Dgroup].fun.Forward(0, x_val);

      auto& iNS = fun[Dgroup].nsRowCol[0];
      auto& matchNS = fun[Dgroup].nsRowCol[2];
      auto& pNS = fun[Dgroup].nsP;

      auto& iFull = fun[Dgroup].outRowCol[0];
      auto& matchFull = fun[Dgroup].outRowCol[2];
      auto& pFull = fun[Dgroup].outP;

//      Rcpp::Rcout << "group " << Dgroup << " " << iNS[0] << "a " << matchNS[0] << "b " << pNS[0] << "c " << iFull[0] << "d " << matchFull[0] << "e " << pFull[0] << "kk" << pFull[1] << "f\n";

      // diagonals
      for(size_t Dk=0; Dk < Nparams; ++Dk) {

        std::fill(direction.begin(), direction.end(), 0.0);
        direction[Dk]  = 1.0;     

        fun[Dgroup].fun.Forward(1, direction);
        fun[Dgroup].fun.Forward(2, directionZeros);

        auto taylor3 = fun[Dgroup].fun.Reverse(3, w);

    // store dense
        const size_t colStart = Dk * Nparams;
        for(size_t Dj=0; Dj<Nparams; Dj++){
          const size_t indexHere = 3*Dj;
          const size_t rowHere = colStart + Dj;
          // first column is T_kii/2 + H_ki
          const double xHere = taylor3[indexHere];//2*taylor3[indexHere]; - Hhere;
          denseThirdHere[rowHere] = xHere;
          if(dense) {
            const double Hhere = taylor3[indexHere+1];
            hessianOutHere[rowHere] += Hhere;
            thirdDiagOutHere[rowHere] += 2*xHere;
          }
        }
      if(!dense) { // not dense
      // fill in the T_kii
      // first column is T_kii/2 + H_ki
        const size_t nsEnd = pNS[Dk+1];
        for(size_t Di=pNS[Dk]; Di<nsEnd; Di++){
          const size_t indexHere = 3*iNS[Di];
          const size_t allHere = matchNS[Di];
          thirdDiagOutHere[allHere] += 2*taylor3[indexHere];
        }
        const size_t fullEnd = pFull[Dk+1];
 
        for(size_t Di=pFull[Dk]; Di<fullEnd; Di++){
          const size_t indexHere = 3*iFull[Di];
          const size_t allHere = matchFull[Di];
          hessianOutHere[allHere] += taylor3[indexHere+1];
        } // Di
      } // else is sparse

      if(Dk == 0) {
        for(size_t Dj=0; Dj<Nparams; Dj++){
          const size_t indexHere = 3*Dj;
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
      TijkHere = std::vector<double>(NtijkHere);

      for(size_t Dpair=0; Dpair < Npairs; ++Dpair) {

        const size_t Di = pairsI[Dpair];
        const size_t Dj = pairsJ[Dpair];
        const size_t pairStart = pairsP[Dpair];
        const size_t DcolHere =  Dpair * Nparams;

//        Rcpp::Rcout << "pair " << Di << " " << Dj << " " << DcolHere << "\n";
        std::fill(direction.begin(), direction.end(), 0.0);
        direction[Di] = direction[Dj] = 1.0;     

        fun[Dgroup].fun.Forward(1, direction);
        fun[Dgroup].fun.Forward(2, directionZeros);
        auto taylor3 = fun[Dgroup].fun.Reverse(3, w);

  // first column is third deriv combination
      //  T_iik/2 + T_jjk/2 + T_ijk 
        const size_t DiNparams = Di* Nparams;
        const size_t DjNparams = Dj* Nparams;
        if(dense) {
          for(size_t Dk=0; Dk<Nparams; Dk++){
            const double TiikTjjk = denseThirdHere[DiNparams + Dk] + denseThirdHere[DjNparams + Dk];
            TijkHere[DcolHere + Dk] = taylor3[3*Dk] - TiikTjjk; 
          }
        } else {
          const size_t NthisPair = pairsPend[Dpair] - pairsP[Dpair];
          for(size_t Dindex=0; Dindex<NthisPair; Dindex++){
            const size_t DindexInIjk = pairStart + Dindex;
            const size_t Dk = sparsityIjk[DindexInIjk];
            const double TiikTjjk = denseThirdHere[DiNparams + Dk] + denseThirdHere[DjNparams + Dk];
            TijkHere[pairStart + Dindex] = taylor3[3*Dk] - TiikTjjk;
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


      for(size_t Dk=0; Dk < Nparams; ++Dk) {
        std::fill(direction.begin(), direction.end(), 0.0);
        direction[Dk]  = 1.0;     
        Qfun.fun.Forward(1, direction);
        Qfun.fun.Forward(2, directionZeros);

        auto taylor3 = Qfun.fun.Reverse(3, w);


      if(dense) {
    // store dense
        const size_t colStart = Dk * Nparams;
        for(size_t Dj=0; Dj<Nparams; Dj++){
          const size_t indexHere = 3*Dj;
          const size_t rowHere = colStart + Dj;

          hessianOutHere[rowHere] += taylor3[indexHere+1];
          thirdDiagOutHere[rowHere] += 2*(taylor3[indexHere]);
        }
      } else { // not dense

        const size_t nsEnd = pNS[Dk+1];
        for(size_t Di=pNS[Dk]; Di<nsEnd; Di++){
          const size_t indexHere = 3*iNS[Di];
          const size_t allHere = matchNS[Di];
          thirdDiagOutHere[allHere] += 2*(taylor3[indexHere]);
        } // Di

        const size_t fullEnd = pFull[Dk+1];
        for(size_t Di=pFull[Dk]; Di<fullEnd; Di++){
          const size_t indexHere = 3*iFull[Di];
          const size_t allHere = matchFull[Di];
          hessianOutHere[allHere] += taylor3[indexHere+1];
        }
      } // else is sparse
      if(Dk == 0) {
        for(size_t Dj=0; Dj<Nparams; Dj++){
          const size_t indexHere = 3*Dj;
          gradientOutHere[Dj] += taylor3[indexHere+2];
        }
      }

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

if (config.verbose ) {
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

/* function for third derivatives, of laplace approx */
//' @export
// [[Rcpp::export]]
Rcpp::List thirdStrata(
  const Rcpp::NumericVector parameters, // beta, gamma, theta
  const Rcpp::List data, 
  const Rcpp::List config,
  SEXP adFun = R_NilValue
  ) {

  const Data   dataC(data);
  const Config configC(config);
  const std::vector<double> x_val(parameters.begin(), parameters.end());

  if(x_val.size() == dataC.Ngamma) Rcpp::warning("need full parameters, only random has been sent");

  AdpackHandle ad = getAdpackFromR(adFun, x_val, dataC, configC);
  std::vector<GroupPack>* fun = ad.ptr;
  auto Qfun = getAdFunQ(x_val, dataC, configC);


//  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];
  const Rcpp::List sparsity = configC.group_sparsity;
  const size_t Ngroup = sparsity.size();

  std::vector<ThirdPack> third(Ngroup);
  if(configC.verbose)
    Rcpp::Rcout << "~";

for (size_t g = 0; g < Ngroup; ++g) {
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
 

//  auto fun = getAdFun(x_val, dataC, configC);

Rcpp::List result = thirdStrata(x_val, third, configC, *fun, Qfun);

return(result);
}
