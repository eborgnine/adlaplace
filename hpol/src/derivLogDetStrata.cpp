#include"hpol.hpp"


struct ThirdPack {
  std::vector<int> pairsI, pairsJ, pairsP, pairsPend, ijkK, ijkMatch;
};


Rcpp::List thirdStrata(
  const CppAD::vector<double>&  parameters, // beta, gamma, theta
  const std::vector<ThirdPack>& third, 
  const Config& config,
  std::vector<GroupPack>& fun,
  GroupPack& Qfun
  ) {


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

  std::vector<double> hessianOut(NoutRowsH, 0.0);
  std::vector<double> thirdDiagOut(NoutRowsT, 0.0);
  std::vector<double> gradientOut(Nparams, 0.0);
  size_t NoutAll = 0;
  for (size_t D = 0; D < Ngroup; ++D) {
    const auto& matchHere = third[D].ijkMatch;
    if (!matchHere.empty()) {
      int m = *std::max_element(matchHere.begin(), matchHere.end());
      if (m >= 0) NoutAll = std::max(NoutAll, static_cast<size_t>(m));
    }
  }
    if (config.verbose ) Rcpp::Rcout << "NoutAll " << NoutAll << "\n";  
  std::vector<double> TijkOutAll(NoutAll+1,0.0);

  if (config.verbose ) Rcpp::Rcout << "starting parallel " << config.num_threads << " threads\n";

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
    CppAD::vector<double> x_val = parameters; 

    std::vector<double> hessianOutHere(NoutRowsH, 0.0);
    std::vector<double> thirdDiagOutHere(NoutRowsT, 0.0);
    std::vector<double> gradientOutHere(Nparams, 0.0);
    // keep all the third Tiik for one group
    std::vector<double> denseThirdHere(NparamsSq, 0.0);

    std::vector<std::vector<double>> TijkThisThread(Ngroup);
    std::vector<size_t> SgroupThisThread(Ngroup);
    size_t DgroupThisThread = 0;

    const std::vector<double> w{0.0, 0.0, 1.0};  
    std::vector<double> direction(Nparams, 0.0);
    const std::vector<double> directionZeros(Nparams, 0.0);
    std::vector<double> taylor3;    


# pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
      SgroupThisThread[DgroupThisThread] = Dgroup;
      std::fill(denseThirdHere.begin(), denseThirdHere.end(), 0.0);

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

        direction[Dk]  = 1.0;     
        fun[Dgroup].fun.Forward(1, direction);
        direction[Dk]  = 0.0;     

        fun[Dgroup].fun.Forward(2, directionZeros);

        std::vector<double>  taylor3tmp = fun[Dgroup].fun.Reverse(3, w);
        taylor3.swap(taylor3tmp);

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
    const size_t NtijkHere = sparsityIjk.size();
    TijkThisThread[DgroupThisThread] = std::vector<double>(NtijkHere);

    for(size_t Dpair=0; Dpair < Npairs; ++Dpair) {

      const size_t Di = pairsI[Dpair];
      const size_t Dj = pairsJ[Dpair];
      const size_t pairStart = pairsP[Dpair];
//      const size_t DcolHere =  Dpair * Nparams;

//        Rcpp::Rcout << "pair " << Di << " " << Dj << " " << DcolHere << "\n";
      direction[Di] = direction[Dj] = 1.0;     
      fun[Dgroup].fun.Forward(1, direction);
      direction[Di] = direction[Dj] = 0.0;     
      fun[Dgroup].fun.Forward(2, directionZeros);
      std::vector<double>  taylor3tmp = fun[Dgroup].fun.Reverse(3, w);
      taylor3.swap(taylor3tmp);


  // first column is third deriv combination
      //  T_iik/2 + T_jjk/2 + T_ijk 
      const size_t DiNparams = Di* Nparams;
      const size_t DjNparams = Dj* Nparams;
      const size_t NthisPair = pairsPend[Dpair] - pairsP[Dpair];
      for(size_t Dindex=0; Dindex<NthisPair; Dindex++){
        const size_t DindexInIjk = pairStart + Dindex;
        const size_t Dk = sparsityIjk[DindexInIjk];
        const double TiikTjjk = denseThirdHere[DiNparams + Dk] + denseThirdHere[DjNparams + Dk];
        TijkThisThread[DgroupThisThread][pairStart + Dindex] =
        taylor3[3*Dk] - TiikTjjk;
      }           
      } // Dpair
      DgroupThisThread++;
  } //group
  const int NgroupThisThread = DgroupThisThread;

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


    // pointer aliases for reduction
  double*       gOut  = gradientOut.data();
  const double* gHere = gradientOutHere.data();
  double*       hOut  = hessianOut.data();
  const double* hHere = hessianOutHere.data();
  double*       tOut  = thirdDiagOut.data();
  const double* tHere = thirdDiagOutHere.data();

    // OpenMP wants an int length
  int nG = static_cast<int>(Nparams);
  int nH = static_cast<int>(NoutRowsH);
  int nT = static_cast<int>(NoutRowsT);

  // scatter into TijkOutAll
    #pragma omp critical(third_accum)
        for (size_t g = 0; g < NgroupThisThread; ++g) {
          const size_t g0   = SgroupThisThread[g];
          const auto& idxes = third[g0].ijkMatch;
          const auto& vals  = TijkThisThread[g];
          for (size_t j = 0; j < idxes.size(); ++j) {
            size_t idx = static_cast<size_t>(idxes[j]);
            TijkOutAll[idx] += vals[j];
          }
  } // g


    #pragma omp critical(grad_accum)
  for (int i = 0; i < nG; ++i) gOut[i] += gHere[i];

    #pragma omp critical(hess_accum)
    for (int i = 0; i < nH; ++i) hOut[i] += hHere[i];

    #pragma omp critical(thirdDiag_accum)
      for (int i = 0; i < nT; ++i) tOut[i] += tHere[i];

} // parallel

if (config.verbose ) {
  Rcpp::Rcout << "done\n";
}

Rcpp::NumericVector TijkR(
  TijkOutAll.begin(), 
  TijkOutAll.end());


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
  const size_t Nparams = parameters.size();
  CppAD::vector<double> x_val(Nparams);

  for(size_t D=0; D<Nparams; D++) {
    x_val[D] = parameters[D];
  }

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
    if (ijk.containsElementNamed("match")) {
      third[g].ijkMatch = Rcpp::as<std::vector<int>>(ijk["match"]);
    } else {
      Rcpp::Rcout << "match missing " << g << "\n";
    }
  }


//  auto fun = getAdFun(x_val, dataC, configC);
  if(configC.verbose)
    Rcpp::Rcout << "-";

  Rcpp::List result = thirdStrata(x_val, third, configC, *fun, Qfun);

  return(result);
}
