#include"hpol.hpp"

#define DEBUG

#define SINGLETHREAD


#include<omp.h>

CppAD::vector<AD<double>> objectiveFunctionInternal(
    CppAD::vector<AD<double>> ad_params, 
    Rcpp::List data, 
    int dirichelet
) {
  

  // S4 objects from data
  Rcpp::S4 QsansDiag(data["QsansDiag"]), A(data["ATp"]),
  X(data["XTp"]), CC(data["cc_matrixTp"]);
  // matrices are transposed, of class dgCMatrix
  
#ifdef DEBUG
    Rcpp::Rcout << "have data, X, A, CC" << std::endl;
#endif
    
  // row, column indices in matrices
  Rcpp::IntegerVector 
  Qrow(QsansDiag.slot("i")), Qcol(QsansDiag.slot("j")),
  Xcol(X.slot("i")), Xp(X.slot("p")), 
  Acol(A.slot("i")), Ap(A.slot("p")), 
  CCcol(CC.slot("i")), CCp(CC.slot("p"));
  // data in matrices
  Rcpp::NumericVector	Qdata =QsansDiag.slot("x"),
    Adata =A.slot("x"), Xdata =X.slot("x");
  
  // number of nonzeros on off-diagonals
  size_t Nq = Qdata.size();
  
  // Q off diagonals
  // separate diagonals and off-diagonals
  Rcpp::NumericVector	Qdiag(data["Qdiag"]);
  
  // map thetas to gammas, y
  Rcpp::IntegerVector map(data["map"]), y(data["y"]);
  
  // dimensions
  Rcpp::IntegerVector dimsX = X.slot("Dim"), dimsA = A.slot("Dim"),
    dimsC = CC.slot("Dim");
  
  // recall X, A, CC transposed
  const size_t Nbeta =dimsX[0], Ngamma = dimsA[0], Neta = dimsA[1],
               Nstrata = dimsC[1];
  const size_t Ntheta = ad_params.size() - Nbeta - Ngamma;
  
  
  CppAD::vector<AD<double>> eta(Neta, 0.0), beta(Nbeta), 
  gamma(Ngamma), theta(Ntheta), logTheta(Ntheta);
  
  
#ifdef DEBUG
    Rcpp::Rcout << "all data" << std::endl;
#endif
    
  // separate out parameters and latent variables
#ifdef DEBUG
  Rcpp::Rcout << "Nbeta " <<  Nbeta << std::endl;
#endif  
  for(size_t D=0;D<Nbeta;D++) {
    beta[D] = ad_params[D];
  }
#ifdef DEBUG
  Rcpp::Rcout << "Ngamma " << Ngamma << std::endl;
#endif
  for(size_t D=0;D<Ngamma;D++) {
    gamma[D] = ad_params[Nbeta+D];
  }
#ifdef DEBUG
  Rcpp::Rcout << "Ntheta " << Ntheta << std::endl;
#endif
  for(size_t D=0;D<Ntheta;D++) {
    theta[D] = ad_params[Nbeta+Ngamma+D];
    logTheta[D] = log(theta[D]);
  }
#ifdef DEBUG
  Rcpp::Rcout << "etas" << std::endl;
#endif
  
  
  // eta = A gamma + X beta
//#pragma omp parallel for
  for(size_t Drow=0; Drow < Neta; Drow++) {
    size_t endCol = Xp[Drow+1];
    for(size_t Dcol=Xp[Drow]; Dcol < endCol; Dcol++) {
      eta[Drow] += Xdata[Dcol] * beta[Xcol[Dcol]];
    }
    endCol = Ap[Drow+1];
    for(size_t Dcol=Ap[Drow]; Dcol < endCol; Dcol++) {
      eta[Drow] += Adata[Dcol] * gamma[Acol[Dcol]];
    }
  }

#ifdef DEBUG
  Rcpp::Rcout << "scaled gammas" << std::endl;
#endif  
  // log(|Q|) + 0.5 * gamma^T Q gamma
  CppAD::vector<AD<double>> gammaScaled(Ngamma);
  AD<double> randomContributionDiag = 0.0; 
  
  // diagonals
  for(size_t D=0;D<Ngamma;D++) {
    size_t mapHere = map[D];
    gammaScaled[D] = gamma[D] / theta[mapHere];
    randomContributionDiag += logTheta[mapHere] +
      0.5*gammaScaled[D]*gammaScaled[D]*Qdiag[D];
  }
  
#ifdef DEBUG
  Rcpp::Rcout << "Q diag" << std::endl;
#endif    
  
  // off diagonals
#ifdef SINGLETHREAD
  std::vector<AD<double>> thread_contrib(1, 0.0);
#else    
  std::vector<AD<double>> thread_contrib(omp_get_max_threads(), 0.0);
#endif  
  
//#pragma omp parallel
{
#ifdef SINGLETHREAD
  size_t thread_num =0; 
#else    
  size_t thread_num = omp_get_thread_num();
#endif  
  AD<double> local_sum = 0.0;
//#pragma omp for
  for(size_t D = 0; D < Nq; D++) {
    local_sum += gammaScaled[Qrow[D]] * gammaScaled[Qcol[D]] * Qdata[D];
  }
  thread_contrib[thread_num] = local_sum;
}

// Combine after parallel region
AD<double> randomContributionOffdiag = 0.0;

#ifdef SINGLETHREAD
randomContributionOffdiag = thread_contrib[0];
#else
  for (auto& val : thread_contrib) randomContributionOffdiag += val;
#endif


#ifdef DEBUG
Rcpp::Rcout << "Q offdiag" << std::endl;
#endif    



// data contribution
AD<double> nu = theta[theta.size()-1];
AD<double> logSqrtNu = log(nu) / 2;
AD<double> oneOverSqrtNu = exp(-logSqrtNu);
AD<double> lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu);

// loop through strata
#ifdef SINGLETHREAD
size_t num_threads = 1;
#else    
size_t num_threads = omp_get_max_threads();
#endif 

std::vector< AD<double> > partial_loglik(num_threads, 0.0);

//CppAD::thread_alloc::parallel_setup(num_threads, in_parallel, thread_number);
//CppAD::thread_alloc::hold_memory(true);
//CppAD::parallel_ad<double>();

//#pragma omp parallel
{
#ifdef SINGLETHREAD
  size_t thread_num = omp_get_thread_num();
#else 
  size_t thread_num =0;
#endif  
  AD<double> thread_loglik = 0.0;
#ifdef DEBUG
  Rcpp::Rcout << "strata loop" << std::endl;
#endif   
//#pragma omp for
  for (size_t i = 0; i < 1; i++){//}Nstrata; i++) {
    
    size_t startHere = CCp[i], Nhere = CCp[i+1];
    AD<double>  contrib = 0.0;
    
    
    // First pass: Compute logSumMu and sumY
    // loop through row i of CCmatrix
    // this is j=startHere
    size_t idx = CCcol[startHere];
    AD<double> logSumMu = eta[idx];
    int sumY = y[idx];
#ifdef DEBUG
    Rcpp::Rcout << startHere << " " << Nhere << " " << logSumMu <<"\n";
#endif 
    for(size_t j=startHere+1; j < Nhere; j++) {
      idx = CCcol[j];
#ifdef DEBUG
      Rcpp::Rcout << " logSumMu " << logSumMu << " idx " <<  idx << 
        " eta[idx] " << eta[idx] << " eta size " << eta.size() <<  " done ";
#endif 
      logSumMu = logspace_add_ad(logSumMu, eta[idx]);
#ifdef DEBUG
      Rcpp::Rcout << "~" << logSumMu << " ";
#endif 
      sumY += y[idx];  		
      }
#ifdef DEBUG
    Rcpp::Rcout << ",";
#endif 
    if(dirichelet) {
      contrib = lgammaOneOverSqrtNu - lgamma_ad(oneOverSqrtNu + sumY);
    } else {
      contrib = 0.0;
    }
#ifdef EVALCONSTANTS
    contrib += lgamma(sumY + 1);
#endif
    
#ifdef DEBUG
    Rcpp::Rcout << ".";
#endif 
    
    // Second pass: Add per-observation terms
    for(size_t j=startHere+1; j < Nhere; j++) {
      size_t idx = CCcol[j];
      AD<double> etaMinusLogSumMu = eta[idx] - logSumMu;
      if(dirichelet) {
        AD<double>  muBarDivSqrtNu = exp(etaMinusLogSumMu - logSqrtNu);
        contrib += lgamma_ad(y[idx] + muBarDivSqrtNu) - lgamma_ad(muBarDivSqrtNu);
      } else {
        contrib += y[idx] * etaMinusLogSumMu;
      }
#ifdef EVALCONSTANTS
      contrib -= lgamma(y[idx] + 1);
#endif
    }
#ifdef DEBUG
    Rcpp::Rcout << ":";
#endif 
    
    thread_loglik += contrib;
  } // loop through strata
#ifdef DEBUG
  Rcpp::Rcout << "strata loop done" << std::endl;
#endif 
  partial_loglik[thread_num] = thread_loglik;
} // parallel region

AD<double> loglik = 0.0;
for (size_t t = 0; t < num_threads; t++) {
  loglik += partial_loglik[t];
}
CppAD::vector<AD<double>> minusLogDens(1,
      -loglik + randomContributionOffdiag + randomContributionDiag);
#ifdef DEBUG
Rcpp::Rcout << "all done" << std::endl;
#endif 
#ifdef EVALCONSTANTS
minusLogDens(0) += Ngamma * HALFLOGTWOPI;
//  minusLogDens -= 0.5*logdet(Q);
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
  
  size_t Nparams = parameters.size();
  int dirichelet = (Rcpp::IntegerVector(config["dirichelet"])[0]);
  CppAD::vector<AD<double>> ad_params(Nparams);
  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  
  if (verbose ) {
    Rcpp::Rcout << "first evaluation:";
  }
  CppAD::vector<AD<double>> y = objectiveFunctionInternal(ad_params, data, dirichelet);
  
  if (verbose ) {
    Rcpp::Rcout << " done." << std::endl;
  }
  CppAD::ADFun<double> f(ad_params, y);
  std::vector<double> x_val(Nparams);
  for(size_t i=0; i < Nparams; i++)
    x_val[i] = parameters[i];
  
  if (verbose ) {
    Rcpp::Rcout << "forward: ";
  }
  
  std::vector<double> y_val = f.Forward(0, x_val);
  
  
  std::vector<double> grad(Nparams);
  grad = f.Jacobian(x_val);
  
  if (verbose ) {
    Rcpp::Rcout << " done." << std::endl;
  }
  
  
  int hesMax; // Default value
  if (config.containsElementNamed("hesMax")) {
    hesMax = config["hesMax"]; // Override if exists
  } else {
    hesMax = parameters.size()*parameters.size()/5; // Default value
  }
  
  int hindex=0;
  Rcpp::NumericVector Hvalue(hesMax);
  Rcpp::IntegerVector Hrow(hesMax), Hcol(hesMax);
  
  std::vector<double> u(Nparams, 0.0);
  std::vector<double> w(1, 1.0);
  std::vector<double> ddw(2*Nparams); 
  double eps = 1e-10, dhere;
  
  if (verbose ) {
    Rcpp::Rcout << "reverse: ";
  }
  
  for (size_t j = 0; j < Nparams; ++j) {
    u[j] = 1.0;
    f.Forward(1, u);
    u[j] = 0.0;  
    
    // Second-order reverse sweep: compute partials of the directional derivative
    ddw = f.Reverse(2, w);
    
    if (verbose ) {
      Rcpp::Rcout << " j";
    }
    for (size_t Drow = 0; Drow <= j; ++Drow) {
      dhere = ddw[2 * Drow + 1];
      if (!CppAD::NearEqual(dhere, 0.0, eps, eps)) {
        if(hindex < hesMax) {
          Hrow[hindex] = Drow;
          Hcol[hindex]= j;
          Hvalue[hindex] = dhere;
        }
        hindex++;
      }
    }
  }
  
  if (verbose ) {
    Rcpp::Rcout << " done." << std::endl;
  }
  
  const int upperBound = std::min(hindex, hesMax);
  
  if(hindex > upperBound)
    Rcpp::warning(
      std::string("increase hesMax, number of hessian elements =  ") + std::to_string(hindex)
    );
  
  
  Rcpp::IntegerVector theDims = Rcpp::IntegerVector::create(Nparams, Nparams);
  
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("value") = y_val,
    Rcpp::Named("grad") = grad,
    Rcpp::Named("hessian") = 
      Rcpp::List::create(
        Rcpp::Named("i") = Rcpp::head(Hrow, upperBound), 
        Rcpp::Named("j") = Rcpp::head(Hcol, upperBound), 
        Rcpp::Named("x") = Rcpp::head(Hvalue, upperBound), 
        Rcpp::Named("dims") = theDims,
        Rcpp::Named("index1") = false,
        Rcpp::Named("symmetric") = true,
        Rcpp::Named("dimnames") = Rcpp::List::create(parameters.names(), parameters.names())
      )
  );
  
  
  return result;
}