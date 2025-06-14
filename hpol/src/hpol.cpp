#include"hpol.hpp"

//#define DEBUG

//#define SINGLETHREAD

#ifndef SINGLETHREAD
#include<omp.h>
#include <vector>
#endif


Rcpp::S4 make_dgTMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N, size_t maxEntries)
{

  int Ni = N;
  Rcpp::IntegerVector dims = 
  Rcpp::IntegerVector::create(Ni, Ni);

    // Convert std::vector<size_t> to IntegerVector (R requires int32 indices)
  Rcpp::NumericVector xR = x[Rcpp::Range(0, maxEntries - 1)];
  Rcpp::IntegerVector iR = i[Rcpp::Range(0, maxEntries - 1)];
  Rcpp::IntegerVector jR = j[Rcpp::Range(0, maxEntries - 1)];


  Rcpp::S4 mat("dgTMatrix");
  mat.slot("i") = iR;
  mat.slot("j") = jR;
  mat.slot("x") = xR;
  mat.slot("Dim") = dims;
//  mat.slot("uplo") =  Rcpp::CharacterVector uplo= Rcpp::wrap('L');
  return mat;
}


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
  
  const std::set<double> unique_map(map.begin(), map.end());
  const size_t Ntheta = unique_map.size() + dirichelet;
  size_t startGamma;


  CppAD::vector<CppAD::AD<double>> beta(Nbeta), 
    gamma(Ngamma), theta(Ntheta), logTheta(Ntheta);
  

#ifdef DEBUG
  Rcpp::Rcout << "Ntheta " << Ntheta << " Neta " << Neta << " Nbeta " << Nbeta <<
  " Ngamma " << Ngamma << "\n";
#endif  

  // eta = A gamma + X beta

  if(ad_params.size() == Ngamma) { // beta and theta are supplied in config
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
  Rcpp::Rcout << "Q offdiag" << std::endl;
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

    for(size_t j=startHere+1; j < Nhere; j++) {
      size_t idx = CCcol[j];
      sumY += y[idx];     
      CppAD::AD<double> etaMinusLogSumMu = eta[idx] - etaLogSum[i];
      CppAD::AD<double>  muBarDivSqrtNu = exp(etaMinusLogSumMu - logSqrtNu);
      if(dirichelet) {
        contrib += lgamma_ad(y[idx] + muBarDivSqrtNu) - lgamma_ad(muBarDivSqrtNu);
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
    Rcpp::Rcout << "eval\n";
  }
  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);
  if (verbose ) {
    Rcpp::Rcout << "eval done " << y[0] << "\n";
  }
  CppAD::ADFun<double> f(ad_params, y);

  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) x_val[i] = parameters[i];

    if (verbose ) {
      Rcpp::Rcout << "forward: ";
    }

    std::vector<double> y_val(1);
    y_val = f.Forward(0, x_val);

    if (verbose ) {
      Rcpp::Rcout << "f0" << y_val[0] << "\n";
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
    std::vector<double> grad = f.Jacobian(x_val);

    result["grad"] = grad;
    if(maxDeriv == 1) {
      return result;
    }

// hessian
    if (verbose ) {
      Rcpp::Rcout << "hess ";
    }

#ifndef SINGLETHREAD
    if (config.containsElementNamed("num_threads")) {
      int num_threads_config = Rcpp::as<int>(config["num_threads"]);
      omp_set_num_threads(num_threads_config);
    }
// Prepare storage for each thread
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif

    std::vector<CppAD::ADFun<double>> f_thread(nthreads);

    for(int D =0; D<nthreads;D++)
      f_thread[D] = f;

//  int hindex=0;
    Rcpp::NumericVector Hvalue(hesMax);
    Rcpp::IntegerVector Hrow(hesMax), Hcol(hesMax);
    double eps = 1e-9;

    std::vector<int> hindex_thread(nthreads, 0);

// Temporary storage for each thread (to avoid race conditions)
    std::vector< std::vector<double> > thread_Hvalue(nthreads, std::vector<double>(hesMax));
    std::vector< std::vector<int> > thread_Hrow(nthreads, std::vector<int>(hesMax));
    std::vector< std::vector<int> > thread_Hcol(nthreads, std::vector<int>(hesMax));

#ifndef SINGLETHREAD
    if (verbose ) {
      Rcpp::Rcout << "parallel\n";
    }


// Parallelized over j (columns)
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      if (verbose ) {
        Rcpp::Rcout <<tid << "\n";
      }
      if (verbose ) {
        Rcpp::Rcout <<tid << "\n";
      }
#else
    if (verbose ) {
      Rcpp::Rcout << "single threaded\n";
    }
      {
        int tid=0;
#endif  

        std::vector<double> u(Nparams, 0.0);
        std::vector<double> w(1, 1.0);
//        std::vector<double> ddw(2*Nparams); 
        double dhere;

        if (verbose ) {
          Rcpp::Rcout <<"thread " << tid << "\n";
        }


        for (int j = tid; j < NparamsI; j+= nthreads) { 
          u[j] = 1.0;    
          auto f1out = f_thread[tid].Forward(1, u);
          u[j] = 0.0;  
          auto ddw = f_thread[tid].Reverse(2, w);

          for (int Drow = 0; Drow < NparamsI; ++Drow) {
            dhere = ddw[2 * Drow + 1];
            if (!CppAD::NearEqual(dhere, 0.0, eps, eps)) {
              if(hindex_thread[tid] < hesMax) {
                thread_Hrow[tid][hindex_thread[tid]] = Drow;
                thread_Hcol[tid][hindex_thread[tid]] = j;
                thread_Hvalue[tid][hindex_thread[tid]] = dhere;
              }
              hindex_thread[tid]++;
            }
          }
        }
        if (verbose ) {
          Rcpp::Rcout <<tid << "\n";
        }
      }

      if (verbose ) {
        Rcpp::Rcout << "combine\n";
      }
// Combine threads' results into the main vectors
      int hindex = 0;
      for (int tid = 0; tid < nthreads; ++tid) {
        for (int k = 0; k < hindex_thread[tid]; ++k) {
          if(hindex < hesMax) {
            Hrow[hindex] = thread_Hrow[tid][k];
            Hcol[hindex] = thread_Hcol[tid][k];
            Hvalue[hindex] = thread_Hvalue[tid][k];
          }
          hindex++;
        }
      }

      if (verbose ) { 
        Rcpp::Rcout << " done." << std::endl;
      }
      const int upperBound = std::min(hindex, hesMax);
      if(hindex > upperBound) 
        Rcpp::warning(
          std::string("increase hesMax, number of hessian elements =  ") + 
          std::to_string(hindex)
          );

      Rcpp::S4 hessianR = make_dgTMatrix(Hvalue, Hrow, Hcol, Nparams, hindex);
      result["hessian"] = hessianR;

      if(verbose) {
        result["hess2"] = Rcpp::List::create(
          Rcpp::Named("x") = Hvalue, Rcpp::Named("i") = Hrow, Rcpp::Named("j") = Hcol,
          Rcpp::Named("nonzeros") = hindex
          );
//      std::vector<double> hess = f.Hessian(x_val, 0);
//      result["hess3"] = hess;
    }

    return result;
}

