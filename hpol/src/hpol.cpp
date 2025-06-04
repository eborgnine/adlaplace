#include"hpol.hpp"

//#define DEBUG


//#define SINGLETHREAD

#ifndef SINGLETHREAD
#include<omp.h>
#endif

#ifdef DEBUG
std::vector<double> out_eta(1000000);
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
    Rcpp::CharacterVector uplo= Rcpp::wrap('L');

    Rcpp::S4 mat("dsTMatrix");
    mat.slot("i") = iR;
    mat.slot("j") = jR;
    mat.slot("x") = xR;
    mat.slot("Dim") = dims;
    mat.slot("uplo") = uplo;
    return mat;
}

atomic_logspace_add& get_logspace_add_atomic_instance() {
    size_t tid = 0;
#ifndef SINGLETHREAD
    tid = omp_get_thread_num();
#endif
    static thread_local atomic_logspace_add instance("logspace_add_atomic_" + std::to_string(tid));
    return instance;
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
    AD<double> logSqrtNu = log(nu) / 2;
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


#ifdef SINGLETHREAD
  size_t max_threads = 1;
#else 
  size_t max_threads = omp_get_max_threads();
#endif  
    
  // off diagonals
  std::vector<AD<double>> thread_offdiagQ(max_threads, 0.0);
  std::vector<AD<double>> thread_loglik(max_threads, 0.0);


#ifndef SINGLETHREAD
#pragma omp parallel
  {
    size_t thread_num = omp_get_thread_num();
    size_t num_threads = omp_get_max_threads();

#else
  {
    size_t thread_num = 0;
    size_t num_threads = 1;
#endif  

//    thread_local atomic_logspace_add logspace_add_atomic_instance("logspace_add_atomic");

//std::string unique_name = "logspace_add_atomic_" + std::to_string(thread_num);

//thread_local atomic_logspace_add logspace_add_atomic_instance(unique_name);

    AD<double> local_offdiagQ = 0.0;
    AD<double> local_loglik = 0.0;


    // Q offdiag    
    for(size_t D = thread_num; D < Nq; D += num_threads) {
        local_offdiagQ += gammaScaled[Qrow[D]] * gammaScaled[Qcol[D]] * Qdata[D];
    }
#ifdef DEBUG
    Rcpp::Rcout << "o " << thread_num << " " <<local_offdiagQ << std::endl;
#endif 

    // loop through strata
    for (size_t i = thread_num; i < Nstrata; 
        i+= num_threads) {

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
        CppAD::vector<AD<double>> in(2), out(1);
        in[0] = logSumMu;
        in[1] = eta[idx];
        get_logspace_add_atomic_instance()(in, out);
        logSumMu = out[0];
// was logSumMu = logspace_add_ad(logSumMu, eta[idx]);

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
      local_loglik += contrib;
  } // loop through strata
  #ifdef DEBUG
    Rcpp::Rcout << "ll " << thread_num  << " " << local_loglik << std::endl;
  #endif  
  thread_offdiagQ[thread_num] = local_offdiagQ;
  thread_loglik[thread_num] = local_loglik;
} // parallel region

#ifdef DEBUG
  Rcpp::Rcout << "strata loop done" << std::endl;
#endif 

// Combine after parallel region
  AD<double> randomContributionOffdiag = 0.0;
  AD<double> loglik = 0.0;
  for (size_t t = 0; t < thread_loglik.size(); t++) {
    randomContributionOffdiag += thread_offdiagQ[t];
    loglik += thread_loglik[t];
  }

#ifdef DEBUG
  Rcpp::Rcout << "logLik " << loglik << " offdiag " << randomContributionOffdiag <<
  " diag " << randomContributionDiag << std::endl;
#endif 

CppAD::vector<AD<double>> minusLogDens(1,
  -loglik + randomContributionOffdiag + 
  randomContributionDiag);

#ifdef EVALCONSTANTS
minusLogDens(0) += Ngamma * HALFLOGTWOPI;
//  minusLogDens -= 0.5*logdet(Q);
#endif
#ifdef DEBUG
Rcpp::Rcout << "all done " << minusLogDens[0] << std::endl;
#endif 

#ifdef DEBUG2
size_t N1 = out_eta.size(), N2 = eta.size();
size_t Nout = (N1 < N2) ? N1 : N2;
for (size_t i = 0; i < Nout; ++i) {
  out_eta[i] = CppAD::Value(eta[i]);
}
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
  
  if (verbose ) {
      Rcpp::Rcout << "eval\n";
  }
  CppAD::vector<AD<double>> y = objectiveFunctionInternal(ad_params, data, dirichelet);
  if (verbose ) {
      Rcpp::Rcout << "eval done " << y[0] << "\n";
  }

  CppAD::ADFun<double> f(ad_params, y);
 
   if (verbose ) {
    Rcpp::Rcout << "1\n";
  }
 std::vector<double> x_val(Nparams);
  if (verbose ) {
    Rcpp::Rcout << "2\n";
  }

  for(size_t i=0; i < Nparams; i++) {
    x_val[i] = parameters[i];
  }
  if (verbose ) {
    Rcpp::Rcout << "forward: ";
  }

  std::vector<double> y_val(1);
  y_val = f.Forward(0, x_val);

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

#ifdef DEBUG2
result["eta"]  = Rcpp::NumericVector(out_eta.begin(), out_eta.end());
#endif

return result;
}

