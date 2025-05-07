#include <TMB.hpp>

template<class Type>
Type loglik_multinom(const vector<Type> &eta,
                     const vector<Type> &y,
                     const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &cc_matrix) {
  Type nll = 0;
  int n_cc = cc_matrix.rows();
  int d_cc = cc_matrix.cols();
  int ccMatrixM1;
  Type lsa, sumY, sumGammaY=0, gammaSumY=0;
  
  for (int i = 0; i < n_cc; i++) {
    lsa = Type(-INFINITY);
    sumY = 0;
    for (int j = 0; j < d_cc; j++) {
      if (cc_matrix(i,j) != 0) {
        ccMatrixM1 = cc_matrix(i,j)-1;
        lsa = logspace_add(lsa, eta(ccMatrixM1));
        sumY +=  y(ccMatrixM1);
        sumGammaY += lgamma(y(ccMatrixM1)+1);
      }
    }
    for (int j = 0; j < d_cc; j++) {
      if (cc_matrix(i,j) != 0) {
        nll += y(cc_matrix(i,j)-1) * (eta(cc_matrix(i,j)-1) - lsa);
      }
    }
    gammaSumY += lgamma(sumY+1);
  }
  
  return nll + gammaSumY - sumGammaY;
}


template<class Type>
Type loglik_dirichlet_multinom(const vector<Type> &eta,
                               const vector<Type> &y,
                               const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &cc_matrix,
                               Type nu) {
  Type nll = 0;
  int n_cc = cc_matrix.rows();
  int d_cc = cc_matrix.cols();
  Type logSumMu, sumY, yHere, logLikConstPart=0, muBarDivSqrtNu;
  Type logSqrtNu = log(nu) / 2;
  Type oneOverSqrtNu = exp(-logSqrtNu);
  Type lgammaOneOverSqrtNu = lgamma(oneOverSqrtNu);
  
  for (int i = 0; i < n_cc; i++) {
    logSumMu = Type(-INFINITY);
    sumY = 0;
    for (int j = 0; j < d_cc; j++) {
      if (cc_matrix(i,j) != 0) {
        logSumMu = logspace_add(logSumMu, eta(cc_matrix(i,j)-1));
        sumY += y(cc_matrix(i,j)-1); 
      }
    }
    
    logLikConstPart += lgammaOneOverSqrtNu + lgamma(1 + sumY) - lgamma(oneOverSqrtNu + sumY);
    
    for (int j = 0; j < d_cc; j++) {
      if (cc_matrix(i,j) != 0) {
        muBarDivSqrtNu = exp(eta(cc_matrix(i,j)-1) - logSumMu - logSqrtNu);
        yHere = y(cc_matrix(i,j)-1);
        nll += lgamma(yHere + muBarDivSqrtNu) - lgamma(yHere + 1) - lgamma(muBarDivSqrtNu);
      }
    }
  }
  
  return nll + logLikConstPart;
}

// Template function
template<class Type>
Type objective_function<Type>::operator() () {
  
  
  // Data inputs
  DATA_SPARSE_MATRIX(X);                   // Design matrix for fixed effects
  DATA_SPARSE_MATRIX(A);                   // Design matrix for random effects
  DATA_VECTOR(y);                   // Observed counts (Poisson responses)
  DATA_SPARSE_MATRIX(Q);            // Full precision matrix
  DATA_IVECTOR(gamma_split);         // Length of gamma for each city
  DATA_VECTOR(psd_scale);         // Length of gamma for each city
  DATA_IMATRIX(cc_matrix);         // Length of gamma for each city
  //  int n_cc = cc_matrix.rows();
  //  int d_cc = cc_matrix.cols();
  
  // Parameter inputs
  PARAMETER_VECTOR(beta);                   // Fixed effects
  PARAMETER_VECTOR(gamma);                  // Concatenated random effects
  PARAMETER_VECTOR(theta);                  // Log-precision parameters (2J values for gamma0 and gamma)
  
  // Compute eta = log(lambda)
  vector<Type> eta(y.size());
  eta.setZero();
  if(X.cols() > 0) eta += X * beta;
  if(A.cols() > 0) eta += A * gamma;
  
  // Compute negative log likelihood (multinomial -- case crossover)
  //  Type lsa;
  
  // components of log density  
  Type logLik, nll = 0, logdetpart=0, randomContribution=0;
  // final parameter is nu
  Type nu = theta(gamma_split.size());
  //  Type logSqrtNu=0, oneOverSqrtNu=0, sumY=0, yHere=0;
  //  Type gammaNu=0, gammaNuSumY=0, muBar=0;
  
  // data contribution
  if (nu > 0)
    nll = loglik_dirichlet_multinom(eta, y, cc_matrix, nu);
  else
    nll = loglik_multinom(eta, y, cc_matrix);
  
  
  // Random effects contribution
  vector<Type> sdTheta(A.cols());
  int idx_start = 0;
  int len_j = 0;
  
  for(int j = 0; j < gamma_split.size(); j++){
    len_j = gamma_split(j);
    // nll -= 0.5 * (len_j*thetaLogPrec(j) + log_dets(j)); // log-determinant part
    //    nll -= 0.5*len_j*theta(j);
    //    nll += len_j*log(theta(j));
    // | sigma^2 Q |^(-1/2) = sigma^P |Q|, logdet = P log sigma + log |Q|
    for(int i = idx_start; i < idx_start + len_j; i++) {
      sdTheta(i) = theta(j)/psd_scale(j); // parameter input is the PSD
      logdetpart += log(sdTheta(i));
    }
    idx_start += len_j;
  }
  
  // do sqrt(exp_theta) * gamma instead of exp_theta * Q
  gamma = gamma * sdTheta;
  randomContribution = -0.5*(gamma * (Q * gamma).col(0)).sum();
  
  logLik = nll + logdetpart + randomContribution;
  REPORT(nll);
  REPORT(logdetpart);
  REPORT(randomContribution);
  REPORT(logLik);
  REPORT(nu);
  
//  Rcout << "; logLik " <<  logLik << " nu " << nu << " nll " << nll;
  
  return -logLik;
}
