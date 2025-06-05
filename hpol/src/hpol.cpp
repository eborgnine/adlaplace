#include"hpol.hpp"

//#define DEBUG

#define SINGLETHREAD

#ifndef SINGLETHREAD
#include<omp.h>
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

  Rcpp::CharacterVector uplo= Rcpp::wrap('L');

  Rcpp::S4 mat("dsTMatrix");
  mat.slot("i") = iR;
  mat.slot("j") = jR;
  mat.slot("x") = xR;
  mat.slot("Dim") = dims;
  mat.slot("uplo") = uplo;
  return mat;
}


CppAD::vector<AD<double>> objectiveFunctionInternal(
  CppAD::vector<AD<double>> ad_params, 
  Rcpp::List data, 
  Rcpp::List config 
  ) {

#ifdef DEBUG
  Rcpp::Rcout << "in objfun\n";
#endif  

  bool verbose = false;
  if (config.containsElementNamed("verbose")) {
    verbose = Rcpp::as<bool>(config["verbose"]);
  }
  bool dirichelet = false;
  if (config.containsElementNamed("dirichelet")) {
    dirichelet = Rcpp::as<bool>(config["dirichelet"]);
  }
  bool transform_theta = false; // theta is logged
  if (config.containsElementNamed("transform_theta")) {
    transform_theta = Rcpp::as<bool>(config["transform_theta"]);
  }

  // S4 objects from data
  Rcpp::S4 QsansDiag(data["QsansDiag"]), A(data["ATp"]),
  X(data["XTp"]), CC(data["cc_matrixTp"]);
  // matrices are transposed, of class dgCMatrix
  
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
  Rcpp::Rcout << "Ntheta " << Ntheta << " Neta " << Neta << " Nbeta " << Nbeta <<
  " Ngamma " << Ngamma << "\n";
#endif  



  // separate out parameters and latent variables
  for(size_t D=0;D<Nbeta;D++) {
    beta[D] = ad_params[D];
  }
  for(size_t D=0;D<Ngamma;D++) {
    gamma[D] = ad_params[Nbeta+D];
  }

  if(transform_theta) {
    for(size_t D=0;D<Ntheta;D++) {
      logTheta[D] = ad_params[Nbeta+Ngamma+D];
      theta[D] = exp(logTheta[D]);      
    }
  } else {
    for(size_t D=0;D<Ntheta;D++) {
      theta[D] = ad_params[Nbeta+Ngamma+D];
      logTheta[D] = log(theta[D]);      
    }
  }

// for data contribution
  AD<double> nu = theta[theta.size()-1];
  AD<double> logSqrtNu = logTheta[theta.size()-1]/ 2;
  AD<double> oneOverSqrtNu = exp(-logSqrtNu);
  AD<double> lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu);


#ifdef DEBUG
  Rcpp::Rcout << "etas";
#endif
  
  
  // eta = A gamma + X beta

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
  Rcpp::Rcout << "." << std::endl;
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
  Rcpp::Rcout << "Q offdiag, strata" << std::endl;
#endif    



  AD<double> local_offdiagQ = 0.0;
  AD<double> local_loglik = 0.0;


    // Q offdiag    
  for(size_t D = 0; D < Nq; D++) {
    local_offdiagQ += gammaScaled[Qrow[D]] * gammaScaled[Qcol[D]] * Qdata[D];
  }

#ifdef DEBUG
  Rcpp::Rcout << "o "  <<local_offdiagQ << std::endl;
#endif 

    // loop through strata
  for (size_t i = 0; i < Nstrata; i++) {

    size_t startHere = CCp[i], Nhere = CCp[i+1];
    AD<double>  contrib = 0.0;

    // First pass: Compute logSumMu and sumY
    // loop through row i of CCmatrix
    // this is j=startHere
    size_t idx = CCcol[startHere];
    AD<double> logSumMu = eta[idx];
    int sumY = y[idx];
    for(size_t j=startHere+1; j < Nhere; j++) {
      idx = CCcol[j];
      sumY += y[idx];     
      logSumMu = logspace_add_ad(logSumMu, eta[idx]);
    }
    if(dirichelet) {
      contrib += lgammaOneOverSqrtNu - lgamma_ad(oneOverSqrtNu + sumY);
    }
#ifdef EVALCONSTANTS
    contrib += lgamma(sumY + 1);
#endif

    // Second pass: Add per-observation terms
    for(size_t j=startHere+1; j < Nhere; j++) {
      size_t idx = CCcol[j];
      AD<double> etaMinusLogSumMu = eta[idx] - logSumMu;
      AD<double>  muBarDivSqrtNu = exp(etaMinusLogSumMu - logSqrtNu);
      if(dirichelet) {
        contrib += lgamma_ad(y[idx] + muBarDivSqrtNu) - lgamma_ad(muBarDivSqrtNu);
      } else {
        contrib += y[idx] * etaMinusLogSumMu;
      }
#ifdef EVALCONSTANTS
      contrib -= lgamma(y[idx] + 1);
#endif
    }
    local_loglik += contrib;
  } // loop through strata
  #ifdef DEBUG
  Rcpp::Rcout << "ll " << local_loglik << std::endl;
  #endif  

#ifdef DEBUG
Rcpp::Rcout << "strata loop done" << std::endl;
#endif 


#ifdef DEBUG
Rcpp::Rcout << "logLik " << local_loglik << " offdiag " << local_offdiagQ <<
" diag " << randomContributionDiag << std::endl;
#endif 

CppAD::vector<AD<double>> minusLogDens(1,
  - local_loglik + local_offdiagQ + 
  randomContributionDiag);

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

  int hesMax = 100000;
  if (config.containsElementNamed("hesMax")) {
    hesMax = Rcpp::as<int>(config["hesMax"]);
  }

  size_t Nparams = parameters.size();


  CppAD::vector<AD<double>> ad_params(Nparams);  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  
  if (verbose ) {
    Rcpp::Rcout << "eval\n";
  }
  CppAD::vector<AD<double>> y = objectiveFunctionInternal(ad_params, data, config);
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

  int hindex=0;
  Rcpp::NumericVector Hvalue(hesMax);
  Rcpp::IntegerVector Hrow(hesMax), Hcol(hesMax);
  
  std::vector<double> u(Nparams, 0.0);
  std::vector<double> w(1, 1.0);
  std::vector<double> ddw(2*Nparams); 
  double dhere, eps = 1e-9;


for (size_t j = 0; j < Nparams; ++j) { 
    u[j] = 1.0;
    f.Forward(1, u);
    u[j] = 0.0;  

    ddw = f.Reverse(2, w);

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
      std::string("increase hesMax, number of hessian elements =  ") + 
      std::to_string(hindex)
    );


  Rcpp::S4 hessianR = make_dgTMatrix(Hvalue, Hrow, Hcol, Nparams, hindex);

  result["hessian"] = hessianR;

#ifdef DEBUG
  result["hes2"] = Rcpp::List::create(
    Rcpp::Named("x") = Hvalue, Rcpp::Named("i") = Hrow, Rcpp::Named("j") = Hcol
    );
#endif  

  return result;
}

