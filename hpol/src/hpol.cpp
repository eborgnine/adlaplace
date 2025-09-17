#include"hpol.hpp"

//#define DEBUG

#include<omp.h>
//#include <vector>


// ---- CppAD parallel hooks ----
bool in_parallel() { return omp_in_parallel() != 0; }
size_t thread_number() { return static_cast<size_t>(omp_get_thread_num()); }



template<class Type>
CppAD::vector<Type>  objectiveFunctionInternal(
  CppAD::vector<Type> ad_params, 
  Rcpp::List data, 
  Rcpp::List config 
  ) {

#ifdef DEBUG
  Rcpp::Rcout << "in objfun\n";
#endif  


  /*  bool verbose = false;
   if (config.containsElementNamed("verbose")) {
   verbose = Rcpp::as<bool>(config["verbose"]);
   }*/

  bool dirichelet = false;
  if (config.containsElementNamed("dirichelet")) {
    dirichelet = Rcpp::as<bool>(config["dirichelet"]);
  }
  bool transform_theta = false; // theta is logged
  if (config.containsElementNamed("transform_theta")) {
    transform_theta = Rcpp::as<bool>(config["transform_theta"]);
  }

  // S4 objects from data
  const Rcpp::S4 QsansDiag(data["QsansDiag"]), A(data["ATp"]),
  X(data["XTp"]), CC(data["cc_matrixTp"]);
  // matrices are transposed, of class dgCMatrix
  
  // row, column indices in matrices
  const Rcpp::IntegerVector 
  Qrow(QsansDiag.slot("i")), Qcol(QsansDiag.slot("j")),
  Xi(X.slot("i")), Xp(X.slot("p")), 
  Ai(A.slot("i")), Ap(A.slot("p")), 
  CCcol(CC.slot("i")), CCp(CC.slot("p"));
  // data in matrices
  const Rcpp::NumericVector	Qdata =QsansDiag.slot("x"),
  Adata =A.slot("x"), Xdata =X.slot("x");
  
  // number of nonzeros on off-diagonals
  const size_t Nq = Qdata.size();
  
  // Q off diagonals
  // separate diagonals and off-diagonals
  const Rcpp::NumericVector	Qdiag(data["Qdiag"]);
  
  // map thetas to gammas, y
  const Rcpp::IntegerVector map(data["map"]), y(data["y"]);
  
  // dimensions
  const Rcpp::IntegerVector dimsX = X.slot("Dim"), dimsA = A.slot("Dim"),
  dimsC = CC.slot("Dim");
  
  // recall X, A, CC transposed
  const size_t Nbeta =dimsX[0], Ngamma = dimsA[0], 
  Neta = dimsA[1], Nstrata = dimsC[1];
//  const size_t NthetaFromParams = ad_params.size() - Nbeta - Ngamma;
  
//  const std::set<int> unique_map(map.begin(), map.end());
  const size_t Ntheta = Rcpp::max(map) + 1 + dirichelet;
  size_t startGamma;

  CppAD::vector<Type> beta(Nbeta), 
  gamma(Ngamma), theta(Ntheta), logTheta(Ntheta);
  

#ifdef DEBUG
  Rcpp::Rcout << "Ntheta " << Ntheta << " Neta " << Neta << " Nbeta " << Nbeta <<
  " Ngamma " << Ngamma << " Dirichelet " << dirichelet << 
  " maxmap " << Rcpp::max(map) << "\n";
#endif  

  // eta = A gamma + X beta

  if(ad_params.size() == Ngamma) { 
  // beta and theta are supplied in config
    startGamma = 0;
    if (config.containsElementNamed("theta")) {
      Rcpp::NumericVector thetaOrig = config["theta"];

      if(transform_theta) {
        for(size_t D=0;D<Ntheta;D++) {
          logTheta[D] = thetaOrig[D];
          theta[D] = exp(logTheta[D]);      
        }
      } else {
        for(size_t D=0;D<Ntheta;D++) {
          theta[D] = thetaOrig[D];
          logTheta[D] = log(theta[D]);      
        }
      }
    } else {
      Rcpp::warning("theta missing from config");
    }

    if (config.containsElementNamed("beta")) {
      Rcpp::NumericVector betaOrig = config["beta"];
      for(size_t D=0; D<Nbeta; D++)
        beta[D] = betaOrig[D];
    } else {
      Rcpp::warning("beta missing from config");    
    }
  } else { // theta and beta in parameters
    if(Ntheta + Ngamma + Nbeta != ad_params.size()) {
      Rcpp::Rcout <<  "parameters is the wrong size: " <<
      "Ntheta " << Ntheta << " Neta " << Neta 
      << " Nbeta " << Nbeta << " Ngamma " << Ngamma << 
      " parameter size " << ad_params.size() << "\n";
    }

    startGamma = Nbeta;
    for(size_t D=0;D<Nbeta;D++) {
      beta[D] = ad_params[D];
    }
    size_t startTheta = Nbeta + Ngamma;
    if(transform_theta) {
      for(size_t D=0;D<Ntheta;D++) {
        logTheta[D] = ad_params[startTheta+D];
        theta[D] = exp(logTheta[D]);      
      }
    } else {
      for(size_t D=0;D<Ntheta;D++) {
        theta[D] = ad_params[startTheta+D];
        logTheta[D] = log(theta[D]);      
      }
    }

  } // end theta and beta in parameters

#ifdef DEBUG
  for(size_t D=0;D<theta.size();D++) {
    Rcpp::Rcout << "D " << D << " theta " << theta[D] << " logtheta " << logTheta[D] << "\n";
  }
#endif          


  // gamma and eta
  for(size_t D=0;D<Ngamma;D++) {
    gamma[D] = ad_params[startGamma+D];
  }

  CppAD::vector<Type> eta(Neta, 0.0);
  
  for(size_t Deta=0; Deta < Neta; Deta++) {
    size_t endX = Xp[Deta+1];
    for(size_t Dbeta=Xp[Deta]; Dbeta < endX; Dbeta++) {
      eta[Deta] += Xdata[Dbeta] * beta[Xi[Dbeta]];
    }
    size_t endA = Ap[Deta+1];
    for(size_t Dgamma=Ap[Deta]; Dgamma < endA; Dgamma++) {
      eta[Deta] += Adata[Dgamma] * gamma[Ai[Dgamma]];
    }
  }


#ifdef DEBUG
  Rcpp::Rcout << "sum log etas\n";
#endif
  // calculate log(sum(exp(eta_i))) within strata
  CppAD::vector<Type> etaLogSum(Nstrata);

  for (size_t i = 0; i < Nstrata; i++) {
    size_t startHere = CCp[i], Nhere = CCp[i+1];
    size_t NinStrata = Nhere - startHere;
    CppAD::vector<Type> etaHere(NinStrata);

    // loop through row i of CCmatrix
    // this is j=startHere
    for(size_t j0=0, j=startHere; j < Nhere; j++,j0++) {
      etaHere[j0] = eta[CCcol[j]];
    }
#ifdef USEATOMICS    
    // use the atomic formula stuff with analytical derivatives
    etaLogSum[i] = logspace_add_n(etaHere);
#else 
    size_t max_idx = 0;
    for(size_t Didx = 1; Didx < etaHere.size(); ++Didx) {
      if(etaHere[Didx] > etaHere[max_idx]) {
        max_idx = Didx;
      }
    }
    double max_value = CppAD::Value(etaHere[max_idx]);
    Type sumexp=0.0;
    for(size_t j = 0; j < etaHere.size(); ++j) {
      sumexp += exp(etaHere[j] - max_value);
    }
    etaLogSum[i] = max_value + log(sumexp);
#endif    
  }




#ifdef DEBUG
  Rcpp::Rcout << "Q diag" << std::endl;
#endif  

  // log(|Q|) + 0.5 * gamma^T Q gamma
  CppAD::vector<Type> gammaScaled(Ngamma);
  Type randomContributionDiag = 0.0; 
  
  // diagonals
  for(size_t D=0;D<Ngamma;D++) {
    size_t mapHere = map[D];

    gammaScaled[D] = gamma[D] / theta[mapHere];
    randomContributionDiag += logTheta[mapHere] +
    0.5*gammaScaled[D]*gammaScaled[D]*Qdiag[D];
  }



#ifdef DEBUG
  Rcpp::Rcout << "Q offdiag " << Nq << std::endl;
#endif    

  Type local_offdiagQ = 0.0;

    // Q offdiag    
  for(size_t D = 0; D < Nq; D++) {
    local_offdiagQ += gammaScaled[Qrow[D]] * gammaScaled[Qcol[D]] * Qdata[D];
  }


#ifdef DEBUG
  Rcpp::Rcout << "data" << std::endl;
#endif    

// for data contribution
  Type nu = theta[theta.size()-1],
  logSqrtNu = logTheta[theta.size()-1]/ 2,
  oneOverSqrtNu = exp(-logSqrtNu),
  lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu),
  local_loglik = 0.0;

#ifdef DEBUG
  Rcpp::Rcout << "nu " << nu << " logSqrtNU " << logSqrtNu << 
  " oneOverSqrtNu " << oneOverSqrtNu << " lgammaOneOverSqrtNu " <<
  lgammaOneOverSqrtNu << std::endl;
#endif    


  // data contribution to loglik, loop through strata
  for (size_t i = 0; i < Nstrata; i++) {

    Type  contrib = 0.0;
    size_t startHere = CCp[i], Nhere = CCp[i+1];
    int sumY = 0;

    for(size_t j=startHere; j < Nhere; j++) {
      size_t idx = CCcol[j];
      sumY += y[idx];     
      Type etaMinusLogSumMu = eta[idx] - etaLogSum[i];
      Type  muBarDivSqrtNu = exp(etaMinusLogSumMu - logSqrtNu);
      if(dirichelet) {
        contrib += lgamma_ad(y[idx] + muBarDivSqrtNu) - 
        lgamma_ad(muBarDivSqrtNu);
      } else {
        contrib += y[idx] * etaMinusLogSumMu;
      }

      
#ifdef EVALCONSTANTS
      contrib -= lgamma(y[idx] + 1);
#endif
    } // j obs in strata

    if(dirichelet) {
      contrib += lgammaOneOverSqrtNu - lgamma_ad(oneOverSqrtNu + sumY);
    }
#ifdef EVALCONSTANTS
    contrib += lgamma(sumY + 1);
#endif

    local_loglik += contrib;
  } // loop through strata


#ifdef DEBUG
  Rcpp::Rcout << "logLik " << local_loglik << " offdiag " << local_offdiagQ <<
  " diag " << randomContributionDiag << std::endl;
#endif 

  CppAD::vector<Type> minusLogDens(1,0);
  minusLogDens[0] =
//  etaLogSum[0]  + etaLogSum[1]
  - local_loglik + local_offdiagQ  + randomContributionDiag
  ;

//    Rcpp::Rcout << " " << etaLogSum[0] << " " << etaLogSum[1] << "\n";

#ifdef EVALCONSTANTS
  minusLogDens[0] += Ngamma * HALFLOGTWOPI;
//  minusLogDens -= 0.5*logdet(Q);
  if (config.containsElementNamed("halfLogDetQ")) {
    minusLogDens[0] -= Rcpp::as<double>(config["halfLogDetQ"]);
  }
#endif

#ifdef DEBUG
  Rcpp::Rcout << "all done " << minusLogDens[0] << std::endl;
#endif 



  return minusLogDens;

}


/* Parameters: 
 X A cc_matrix: dgRMatrix, cc_matrix can be lgRMatrix or ngRMatrix (with MatrixExtra package)
 QsansDiag, Qdiag dsTMatrix and vector for precision matrix
 map, y integer vectors 
 config: hesMax integer, maximum number of non-zero hessian elements
 */
//' @export
// [[Rcpp::export]]
Rcpp::List objectiveFunctionC(
  Rcpp::NumericVector parameters, 
  Rcpp::List data, 
  Rcpp::List config
  ) {

  bool verbose = false;
  if (config.containsElementNamed("verbose")) {
    verbose = Rcpp::as<bool>(config["verbose"]);
  }
  size_t maxDeriv = 2;
  if (config.containsElementNamed("maxDeriv")) {
    maxDeriv = Rcpp::as<int>(config["maxDeriv"]);
  }
  int num_threads = 1;
  if (config.containsElementNamed("num_threads"))
    num_threads = Rcpp::as<int>(config["num_threads"]);


  size_t Nparams = parameters.size();
  int NparamsI = static_cast<int>(Nparams);
  int hesMax = NparamsI*NparamsI/2;
  if (config.containsElementNamed("hesMax")) {
    hesMax = Rcpp::as<int>(config["hesMax"]);
  }



  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  
  if (verbose ) {
    Rcpp::Rcout << "eval";
  }
  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);
  if (verbose ) {
    Rcpp::Rcout << "y " << y[0] << "\n";
  }
  CppAD::ADFun<double> fun(ad_params, y);

  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) { 
    x_val[i] = parameters[i];
  }

  if (verbose ) {
    Rcpp::Rcout << "forward: ";
  }

  std::vector<double> y_val(1);
  y_val = fun.Forward(0, x_val);

  if (verbose ) {
    Rcpp::Rcout << "f0 " << y_val[0] << "\n";
  }


  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("value") = y_val[0]
    );

  if(maxDeriv == 0) {
    return result;
  }

  if (verbose ) {
    Rcpp::Rcout << "grad ";
  }
  std::vector<double> grad = fun.Jacobian(x_val);
  if (verbose ) {
    Rcpp::Rcout << ".\n";
  }
  result["grad"] = grad;
  if(maxDeriv == 1) {
    return result;
  }

// hessian
  if (verbose ) {
    Rcpp::Rcout << "hess " << num_threads << " threads\n";
  }

  CppAD::vector< std::set<size_t> >  sparsity(Nparams);
  Rcpp::IntegerVector Hrow, Hcol, Hp;
  Rcpp::NumericVector Hvalue;
  bool haveSparsity;

  Rcpp::IntegerVector Hstart(num_threads,0), Hend(num_threads,0);


#ifdef DEBUG
  Rcpp::NumericMatrix denseHessianOut(Nparams, Nparams);
#endif

    // if have sparsity
  if (config.containsElementNamed("sparsity")) {
    if (verbose ) {
      Rcpp::Rcout << "given sparsity pattern\n";
    }      
    haveSparsity=TRUE;

    Rcpp::List sparsityR, sparsityR2, sparsityR1 = config["sparsity"];
    if (sparsityR1.containsElementNamed("second")) {
        // list with second and third sparsity
      sparsityR2 = sparsityR1["second"];
      if (config.containsElementNamed("theta")) {
        if (verbose ) {
          Rcpp::Rcout << "sparsity only for random effects\n";
        }   
        sparsityR = sparsityR2["random"];
      } else {
        if (verbose ) {
          Rcpp::Rcout << "full parameter sparsity\n";
        }   
        sparsityR = sparsityR2["full"];
      }
    } else {
      if (verbose ) {
        Rcpp::Rcout << "sparsity i j p provided\n";
      }   
      sparsityR = sparsityR1;
    }

    Hrow = sparsityR["i"]; 
    Hcol = sparsityR["j"];
    Hp = sparsityR["p"];
    Hvalue = Rcpp::NumericVector(Hrow.size());


    // where each thread starts and ends in hessian
    int NperThread = (Hvalue.size() + num_threads - 1) / num_threads;

    if (verbose ) {
      Rcpp::Rcout << "NperThread " << NperThread << 
      " non-zeros " << Hvalue.size() << "\n";
    }   


  // start and end points for each thread
    Hstart = NperThread * Rcpp::seq(0, num_threads - 1);
  Hstart = Rcpp::pmin(Hstart, Hrow.size());  // clamp
  Hend = Rcpp::pmin(Hstart + NperThread, Hrow.size());
  



      // convert to vector set
// Loop over columns (variables)
  for(int D=0;D<Hrow.size();++D) {
    int Drow = Hrow[D], Dcol=Hcol[D];
    sparsity[Drow].insert(Dcol);    
    if(Drow != Dcol) sparsity[Dcol].insert(Drow);   
  }


    } else { // if don't have sparsity
      if (verbose ) {
        Rcpp::Rcout << "compute dense hessian hesMax " << hesMax << "\n";
      }
      haveSparsity = FALSE;
      Hrow = Rcpp::IntegerVector(hesMax, -1); 
      Hcol = Rcpp::IntegerVector(hesMax, -1); 
      Hvalue = Rcpp::NumericVector(hesMax);
      int NperThread = Hvalue.size()/num_threads;
      Hend = Hstart = NperThread * Rcpp::seq(0, num_threads - 1); 
    }



    omp_set_num_threads(num_threads);
    CppAD::thread_alloc::parallel_setup(num_threads, in_parallel, thread_number);
    CppAD::thread_alloc::hold_memory(true);

    // Replicate fun object for each thread
    std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
    for (int i = 0; i < num_threads; ++i) {
      fun_threads[i] = fun;
    }

    const double eps = 1e-9;
//    int hindex = 0, hindex_thread=0;



    #pragma omp parallel
    {


      const int tid = thread_number();
      const int nthreads_thread = omp_get_num_threads();
      std::vector<double> w(1, 1.0);

      // 
      if(haveSparsity) {
        // rounding up
        int Nelements = Hend[tid] - Hstart[tid];
        std::vector<double> thread_Hvalue(Nelements);
        std::vector<int> thread_Hrow(Nelements), thread_Hcol(Nelements);

        for(int D=0,Dorig=Hstart[tid];D<Nelements;D++,Dorig++) {
          thread_Hrow[D] = Hrow[Dorig];
          thread_Hcol[D] = Hcol[Dorig];
        }

#ifdef DEBUG
        Rcpp::Rcout << "t" << tid << "N" << Nelements <<
        "s" << Hstart[tid] << "e" << Hend[tid] <<
        "\n";
#endif 

        CppAD::sparse_hessian_work work; 

        fun_threads[tid].SparseHessian(
          x_val, w, sparsity, 
          thread_Hrow, thread_Hcol, thread_Hvalue, 
          work);

// copy elements to result vector
        for(
          int DfromZero=0,DfromStart=Hstart[tid];
          DfromZero < thread_Hvalue.size(); //DfromStart < Hend[tid]; 
          ++DfromZero,++DfromStart
          ) {
          Hvalue[DfromStart] = thread_Hvalue[DfromZero];
#ifdef DEBUG
        denseHessianOut(thread_Hrow[DfromZero], 
          thread_Hcol[DfromZero]) = thread_Hvalue[DfromZero];
#endif

      }

      } else { // no sparsity
        // eval dense hessian, save non-zeros


#ifdef DEBUG
        Rcpp::Rcout << "t" << tid << "\n";
#endif 

        std::vector<double> u(Nparams, 0.0);

    #pragma omp for
        for (int j = 0; j < NparamsI; j++) {
          std::fill(u.begin(), u.end(), 0.0);
          u[j] = 1.0;
          fun_threads[tid].Forward(0, x_val);
          fun_threads[tid].Forward(1, u);
          std::vector<double> ddw = fun_threads[tid].Reverse(2, w);

          for (int irow = j; irow < NparamsI; ++irow) {
            double dhere = ddw[2 * irow + 1];
#ifdef DEBUG
            denseHessianOut(irow, j) = dhere;
#endif
            if (!CppAD::NearEqual(dhere, 0.0, eps, eps)) {
              Hvalue[Hend[tid]] = dhere;
              Hrow[Hend[tid]] = irow;
              Hcol[Hend[tid]] = j;
              if(Hend[tid]<hesMax) Hend[tid]++;
            } // if dhere not zero
          } // irow
      } // j column
    } // no sparsity pattern

  } // parallel


#ifdef UNDEF
  CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
  CppAD::thread_alloc::hold_memory(false);
  for (int i = 1; i < num_threads; ++i)
    CppAD::thread_alloc::free_available(i);
  CppAD::thread_alloc::free_available(0);
#endif

/*Rcpp::List hessianR = Rcpp::List::create(
  Rcpp::Named("i") = Hrow,
  Rcpp::Named("j") = Hcol,
  Rcpp::Named("x") = Hvalue,
  Rcpp::Named("start") = Hstart,
  Rcpp::Named("end") = Hend
);*/

  Rcpp::S4 hessianR = make_TMatrix(
    Hvalue, Hrow, Hcol, 
    Nparams); 

  result["hessian"] = hessianR;

#ifdef DEBUG
  result["denseHessian"] = denseHessianOut;
#endif

  if (verbose ) { 
    Rcpp::Rcout << " done." << std::endl;
  }

  return result;
}

