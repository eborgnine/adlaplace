#include"hpol.hpp"

#define DEBUG

#include<omp.h>
//#include <vector>


// ---- CppAD parallel hooks ----
bool in_parallel() { return omp_in_parallel() != 0; }
size_t thread_number() { return static_cast<size_t>(omp_get_thread_num()); }



CppAD::vector<CppAD::AD<double>> objectiveFunctionInternal(
  CppAD::vector<CppAD::AD<double>> ad_params, 
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
  
  const std::set<int> unique_map(map.begin(), map.end());
  const size_t Ntheta = unique_map.size() + dirichelet;
  size_t startGamma;

  CppAD::vector<CppAD::AD<double>> beta(Nbeta), 
  gamma(Ngamma), theta(Ntheta), logTheta(Ntheta);
  

#ifdef DEBUG
  Rcpp::Rcout << "Ntheta " << Ntheta << " Neta " << Neta << " Nbeta " << Nbeta <<
  " Ngamma " << Ngamma << "\n";
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
      Rcpp::warning("parameters is the wrong size");
      Rcpp::Rcout << "Ntheta " << Ntheta << " Neta " << Neta 
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


  // gamma and eta
  for(size_t D=0;D<Ngamma;D++) {
    gamma[D] = ad_params[startGamma+D];
  }

  CppAD::vector<CppAD::AD<double>> eta(Neta, 0.0);
  
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
  CppAD::vector<CppAD::AD<double>> etaLogSum(Nstrata);

  for (size_t i = 0; i < Nstrata; i++) {
    size_t startHere = CCp[i], Nhere = CCp[i+1];
    size_t NinStrata = Nhere - startHere;
    CppAD::vector<CppAD::AD<double>> etaHere(NinStrata);

    // loop through row i of CCmatrix
    // this is j=startHere
    for(size_t j0=0, j=startHere; j < Nhere; j++,j0++) {
      etaHere[j0] = eta[CCcol[j]];
    }
    etaLogSum[i] = logspace_add_n(etaHere);
  }




#ifdef DEBUG
  Rcpp::Rcout << "Q diag" << std::endl;
#endif  

  // log(|Q|) + 0.5 * gamma^T Q gamma
  CppAD::vector<CppAD::AD<double>> gammaScaled(Ngamma);
  CppAD::AD<double> randomContributionDiag = 0.0; 
  
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

  CppAD::AD<double> local_offdiagQ = 0.0;

    // Q offdiag    
  for(size_t D = 0; D < Nq; D++) {
    local_offdiagQ += gammaScaled[Qrow[D]] * gammaScaled[Qcol[D]] * Qdata[D];
  }


#ifdef DEBUG
  Rcpp::Rcout << "data" << std::endl;
#endif    

// for data contribution
  CppAD::AD<double> nu = theta[theta.size()-1];
  CppAD::AD<double> logSqrtNu = logTheta[theta.size()-1]/ 2;
  CppAD::AD<double> oneOverSqrtNu = exp(-logSqrtNu);
  CppAD::AD<double> lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu);




  // data contribution to loglik, loop through strata
  CppAD::AD<double> local_loglik = 0.0;
  for (size_t i = 0; i < Nstrata; i++) {

    CppAD::AD<double>  contrib = 0.0;
    size_t startHere = CCp[i], Nhere = CCp[i+1];
    int sumY = 0;

    for(size_t j=startHere; j < Nhere; j++) {
      size_t idx = CCcol[j];
      sumY += y[idx];     
      CppAD::AD<double> etaMinusLogSumMu = eta[idx] - etaLogSum[i];
      CppAD::AD<double>  muBarDivSqrtNu = exp(etaMinusLogSumMu - logSqrtNu);
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

  CppAD::vector<CppAD::AD<double>> minusLogDens(1,0);
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
  int hesMax = NparamsI*NparamsI;
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
    Rcpp::Rcout << " " << y[0] << "\n";
  }
  CppAD::ADFun<double> fun(ad_params, y);

  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) x_val[i] = parameters[i];

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

    Rcpp::IntegerVector Hrow, Hcol;
    bool haveSparsity;
    CppAD::vector< std::set<size_t> >  sparsity(Nparams);

    if (config.containsElementNamed("sparsity")) {
      if (verbose ) {
        Rcpp::Rcout << "given sparsity pattern\n";
      }      
      haveSparsity=TRUE;
      Rcpp::List sparsityR = config["sparsity"];
      Hrow = sparsityR["i"];
      Hcol = sparsityR["j"];
      auto Hp = compute_p_vector(Hcol, Nparams);//sparsityR["p"];
// Loop over columns (variables)
      for (size_t col = 0; col < Nparams; ++col) {
        for (int idx = Hp[col]; idx < Hp[col + 1]; ++idx) {
          size_t row = Hrow[idx];
          sparsity[col].insert(row);
        // For symmetric sparsity (Hessian), also do:
          sparsity[row].insert(col);
        }
      }
    } else {
#ifdef DEBUG
      Rcpp::Rcout << "no sparsity pattern\n";
#endif      
      Hrow = Rcpp::IntegerVector(hesMax);
      Hcol = Rcpp::IntegerVector(hesMax);
      haveSparsity=FALSE;
    }

    Rcpp::NumericVector Hvalue(Hrow.size());

    
    omp_set_num_threads(num_threads);
    CppAD::thread_alloc::parallel_setup(num_threads, in_parallel, thread_number);
    CppAD::thread_alloc::hold_memory(true);

    // Replicate fun object for each thread
    std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
    for (int i = 0; i < num_threads; ++i) fun_threads[i] = fun;

      const double eps = 1e-9;
    int hindex = 0;

// TO DO: https://cppad.readthedocs.io/latest/sparse_hessian.html
    // https://cppad.readthedocs.io/latest/RevSparseHes.html
    // keep sparsity pattern, get pairs of row, col, divide into equal parts

    #pragma omp parallel
    {    


      const int tid = thread_number();
      const int nthreads_thread = omp_get_num_threads();
      int hindex_thread = 0, Hstart=0;
      std::vector<double> w(1, 1.0);
      std::vector<double> thread_Hvalue;
      std::vector<int> thread_Hrow, thread_Hcol;
        int NperThread = Hvalue.size() / nthreads_thread;
        if(NperThread * nthreads_thread < Hvalue.size()) {
          NperThread += 1;
        }

      // 
      if(haveSparsity){
        Hstart = NperThread * tid;
        int Hend = Hstart + NperThread;
        if(Hend > Hvalue.size()) Hend = Hvalue.size();
        hindex_thread = Hend - Hstart;



        CppAD::sparse_hessian_work work; 

        thread_Hrow = std::vector<int>(
          Hrow.begin() + Hstart, Hrow.begin() + Hend);
        thread_Hcol = std::vector<int>(
          Hcol.begin() + Hstart, Hcol.begin() + Hend);
        thread_Hvalue.resize(hindex_thread);
   
#ifdef DEBUG
      Rcpp::Rcout << "t" << tid << "n"<< thread_Hvalue.size() << "\n";
#endif 

        fun_threads[tid].SparseHessian(
          x_val, w, sparsity, 
          thread_Hrow, thread_Hcol, thread_Hvalue, 
          work);

      } else {
        // eval dense hessian, save non-zeros
        thread_Hrow.resize(NperThread);
        thread_Hcol.resize(NperThread);
        thread_Hvalue.resize(NperThread);
        std::vector<double> u(Nparams, 0.0);

        for (int j = tid; j < NparamsI; j += nthreads_thread) {
          std::fill(u.begin(), u.end(), 0.0);
          u[j] = 1.0;
          fun_threads[tid].Forward(0, x_val);
          fun_threads[tid].Forward(1, u);
          std::vector<double> ddw = fun_threads[tid].Reverse(2, w);
          for (int irow = 0; irow < NparamsI; ++irow) {
            double dhere = ddw[2 * irow + 1];
            if (!CppAD::NearEqual(dhere, 0.0, eps, eps)) {
              if (hindex_thread < thread_Hrow.size()) {
                thread_Hrow[hindex_thread] = irow;
                thread_Hcol[hindex_thread] = j;
                thread_Hvalue[hindex_thread] = dhere;
                hindex_thread++;
              }
            }
          }
      } // j column
    } // no sparsity pattern

    #pragma omp critical
    {
      int out_start = haveSparsity ? Hstart : hindex;
      for (int k = 0; k < hindex_thread; ++k) {
        if(!haveSparsity) {
          Hrow[out_start + k] = thread_Hrow[k];
          Hcol[out_start + k] = thread_Hcol[k];
        }
        Hvalue[out_start + k] = thread_Hvalue[k];
      }
      hindex += hindex_thread;
    }
  }
  if (verbose ) { 
    Rcpp::Rcout << "hessian " << hindex << " entries" << std::endl;
  }

  CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
  CppAD::thread_alloc::hold_memory(false);
  for (int i = 1; i < num_threads; ++i)
    CppAD::thread_alloc::free_available(i);
  CppAD::thread_alloc::free_available(0);

  Rcpp::S4 hessianR = make_dgTMatrix(Hvalue, Hrow, Hcol, Nparams, hindex);
  result["hessian"] = hessianR;


  if (verbose ) { 
    Rcpp::Rcout << " done." << std::endl;
  }

  if(verbose) {
    result["hess2"] = Rcpp::List::create(
      Rcpp::Named("x") = Hvalue, Rcpp::Named("i") = Hrow, Rcpp::Named("j") = Hcol,
      Rcpp::Named("nonzeros") = hindex
      );

  }

  return result;
}

