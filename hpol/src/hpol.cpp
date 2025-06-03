#include"hpol.hpp"

//#define DEBUG

#define SINGLETHREAD

#ifdef SINGLETHREAD
#include<omp.h>
#endif

Rcpp::S4 make_dgTMatrix(
    const std::vector<double>& x,
    const std::vector<size_t>& i,
    const std::vector<size_t>& j,
    int N)
{

    Rcpp::IntegerVector dims = 
      Rcpp::IntegerVector::create(N, N);

    // Convert std::vector<size_t> to IntegerVector (R requires int32 indices)
    Rcpp::IntegerVector iR(i.begin(), i.end());
    Rcpp::IntegerVector jR(j.begin(), j.end());
    Rcpp::NumericVector xR(x.begin(), x.end());

    Rcpp::S4 mat("dgTMatrix");
    mat.slot("i") = iR;
    mat.slot("j") = jR;
    mat.slot("x") = xR;
    mat.slot("Dim") = dims;
    return mat;
}

CppAD::vector<AD<double>> objectiveFunctionInternal(
  CppAD::vector<AD<double>> ad_params, 
  Rcpp::List data, 
  int dirichelet
  ) {


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
  
#ifdef SINGLETHREAD
  size_t num_threads = 1;
#else    
  size_t num_threads = omp_get_max_threads();
#endif 

  // separate out parameters and latent variables
  for(size_t D=0;D<Nbeta;D++) {
    beta[D] = ad_params[D];
  }
  for(size_t D=0;D<Ngamma;D++) {
    gamma[D] = ad_params[Nbeta+D];
  }

  for(size_t D=0;D<Ntheta;D++) {
    theta[D] = ad_params[Nbeta+Ngamma+D];
    logTheta[D] = log(theta[D]);
  }

#ifdef DEBUG
  Rcpp::Rcout << "etas" << std::endl;
#endif
  
  
  // eta = A gamma + X beta
#ifndef SINGLETHREAD  
#pragma omp parallel for
#endif
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
  
  
  // off diagonals
  std::vector<AD<double>> thread_contrib(num_threads, 0.0);
  
#ifndef SINGLETHREAD
#pragma omp parallel
  {
    size_t thread_num = omp_get_thread_num();
#else
    {
      size_t thread_num = 0;
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
    for (size_t t = 0; t < num_threads; t++) {
      randomContributionOffdiag += thread_contrib[t];
    }


#ifdef DEBUG
    Rcpp::Rcout << "Q offdiag" << std::endl;
#endif    



// data contribution
    AD<double> nu = theta[theta.size()-1];
    AD<double> logSqrtNu = log(nu) / 2;
    AD<double> oneOverSqrtNu = exp(-logSqrtNu);
    AD<double> lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu);

// loop through strata

    std::vector< AD<double> > partial_loglik(num_threads, 0.0);


#ifndef SINGLETHREAD
#pragma omp parallel
    {
      size_t thread_num = omp_get_thread_num();
# else 
      {    
        size_t thread_num = 0;
#endif  
        AD<double> thread_loglik = 0.0;
#ifdef DEBUG
        Rcpp::Rcout << "strata loop" << std::endl;
#endif  
#ifndef SINGLETHREAD 
#pragma omp for
#endif
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
            CppAD::vector<AD<double>> in(2), out(1);
            in[0] = logSumMu;
            in[1] = eta[idx];
            logspace_add_atomic(in, out);
            logSumMu = out[0];
//        logSumMu = logspace_add_ad(logSumMu, eta[idx]);
            sumY += y[idx];  		
          }
          if(dirichelet) {
            contrib = lgammaOneOverSqrtNu - lgamma_ad(oneOverSqrtNu + sumY);
          } else {
            contrib = 0.0;
          }
#ifdef EVALCONSTANTS
          contrib += lgamma(sumY + 1);
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
  -loglik + randomContributionOffdiag + 
  randomContributionDiag);

#ifdef EVALCONSTANTS
minusLogDens(0) += Ngamma * HALFLOGTWOPI;
//  minusLogDens -= 0.5*logdet(Q);
#endif
#ifdef DEBUG
Rcpp::Rcout << "all done" << std::endl;
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
  int dirichelet = (Rcpp::IntegerVector(config["dirichelet"])[0]);
  CppAD::vector<AD<double>> ad_params(Nparams);
  
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  

  CppAD::vector<AD<double>> y = objectiveFunctionInternal(ad_params, data, dirichelet);
  

    CppAD::ADFun<double> f(ad_params, y);
    std::vector<double> x_val(Nparams);
    for(size_t i=0; i < Nparams; i++)
      x_val[i] = parameters[i];

    if (verbose ) {
      Rcpp::Rcout << "forward: ";
    }

    std::vector<double> y_val = f.Forward(0, x_val);

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

  // Before RevSparseHes, you MUST call ForSparseJac with a set-based pattern
  CppAD::vector<std::set<size_t>> jac_pattern(Nparams);
  for (size_t i = 0; i < Nparams; ++i)
    jac_pattern[i].insert(i); // Identity pattern

  f.ForSparseJac(Nparams, jac_pattern);  // <-- REQUIRED!


  CppAD::vector< std::set<size_t> > select_range(1);//, hes_pattern(Nparams);

  select_range[0].insert(0);
  bool transpose = false;

auto hes_pattern = f.RevSparseHes(Nparams, select_range, transpose);



std::vector<size_t> row, col;
for (size_t i = 0; i < Nparams; ++i) {
     std::set<size_t>::const_iterator itr;
     for(itr = hes_pattern[i].begin(); itr != hes_pattern[i].end(); itr++) {
          auto j = *itr;
          if (j <= i) { 
            row.push_back(i);
            col.push_back(j);
        }
      }
    }

  if (verbose ) {
    Rcpp::Rcout << ".";
  }
CppAD::vector<double> cppad_x(x_val.size()), cppad_w(1, 1.0), hes(row.size());
CppAD::vector<size_t> cppad_row(row.size()), cppad_col(col.size());
CppAD::vector<std::set<size_t>> cppad_pattern(hes_pattern.size());

for (size_t i = 0; i < x_val.size(); ++i) cppad_x[i] = x_val[i];
for (size_t i = 0; i < row.size(); ++i) cppad_row[i] = row[i];
for (size_t i = 0; i < col.size(); ++i) cppad_col[i] = col[i];
for (size_t i = 0; i < hes_pattern.size(); ++i) cppad_pattern[i] = hes_pattern[i];
  if (verbose ) {
    Rcpp::Rcout << ".";
  }
CppAD::sparse_hessian_work work;
f.SparseHessian(cppad_x, cppad_w, cppad_pattern, cppad_row, cppad_col, hes, work);

  if (verbose ) {
    Rcpp::Rcout << ".";
  }
// Copy to STL/Rcpp for return
std::vector<double> hes_x(hes.size());
std::copy(hes.begin(), hes.end(), hes_x.begin());


Rcpp::S4 hessianR = make_dgTMatrix(hes_x, row, col, Nparams);

result["hessian"] = hessianR;

return result;
}

