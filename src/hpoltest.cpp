#include <TMB.hpp>
#define EIGEN_DONT_PARALLELIZE

//#define EVALCONSTANTS
#define LOGTWOPI 1.8378770664093454835606594728112352797227949472755668


template<class Type>
Type loglik_multinom(const vector<Type> &eta,
                     const vector<Type> &y,
                     const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &cc_matrix) {
  Type loglik = 0;
  int n_cc = cc_matrix.rows();
  int d_cc = cc_matrix.cols();


#ifdef EVALCONSTANTS
  Type sumGammaY=0, gammaSumY=0, sumY;
#endif
  for (int i = 0; i < n_cc; i++) PARALLEL_REGION {
    Type lsa = Type(-INFINITY);
    int ccMatrixM1;
#ifdef EVALCONSTANTS   
    Type sumY = 0;
#endif
    for (int j = 0; j < d_cc; j++) {
      ccMatrixM1 = cc_matrix(i,j);
      if (ccMatrixM1 != 0) {
        ccMatrixM1 -= 1;
        lsa = logspace_add(lsa, eta(ccMatrixM1));
#ifdef EVALCONSTANTS        
        sumY +=  y(ccMatrixM1);
        sumGammaY += lgamma(y(ccMatrixM1)+1);
#endif        
      }
    }
    for (int j = 0; j < d_cc; j++) {
      ccMatrixM1 = cc_matrix(i,j);
      if (ccMatrixM1 != 0) {
        ccMatrixM1 -= 1;
        loglik += y(ccMatrixM1) * (eta(ccMatrixM1) - lsa);
      }
    }
#ifdef EVALCONSTANTS        
    gammaSumY += lgamma(sumY+1);
#endif
    }
#ifdef EVALCONSTANTS  
  loglik += gammaSumY - sumGammaY;
#endif
  return loglik;
}


template<class Type>
Type loglik_dirichlet_multinom(const vector<Type> &eta,
                               const vector<Type> &y,
                               const Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> &cc_matrix,
                               Type nu) {
  Type loglik = 0;
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
    
    logLikConstPart += lgammaOneOverSqrtNu  - lgamma(oneOverSqrtNu + sumY);
#ifdef EVALCONSTANTS    
    logLikConstPart += lgamma(1 + sumY);
#endif    
    for (int j = 0; j < d_cc; j++) {
      if (cc_matrix(i,j) != 0) {
        muBarDivSqrtNu = exp(eta(cc_matrix(i,j)-1) - logSumMu - logSqrtNu);
        yHere = y(cc_matrix(i,j)-1);
        loglik += lgamma(yHere + muBarDivSqrtNu)  - lgamma(muBarDivSqrtNu);
#ifdef EVALCONSTANTS        
        loglik -= lgamma(yHere + 1);
#endif          
      }
    }
  }
  
  return loglik + logLikConstPart;
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
  Type loglik = 0, logdetpart=0, randomContribution=0;
  // final parameter is nu
  Type nu = theta(gamma_split.size());
  //  Type logSqrtNu=0, oneOverSqrtNu=0, sumY=0, yHere=0;
  //  Type gammaNu=0, gammaNuSumY=0, muBar=0;
  
  // data contribution
  if (nu > 0) {
    loglik = loglik_dirichlet_multinom(eta, y, cc_matrix, nu);
  } else {
    loglik = loglik_multinom(eta, y, cc_matrix);
  }
  
  // Random effects contribution
  vector<Type> sdTheta(A.cols());
  int idx_start = 0;
  int len_j = 0;
  
  for(int j = 0; j < gamma_split.size(); j++){
    len_j = gamma_split(j);
    // loglik -= 0.5 * (len_j*thetaLogPrec(j) + log_dets(j)); // log-determinant part
    //    loglik -= 0.5*len_j*theta(j);
    //    loglik += len_j*log(theta(j));
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
#ifdef EVALCONSTANTS
  randomContribution += (gamma.size()/2)*LOGTWOPI
#endif  
  
  loglik += logdetpart + randomContribution;
  REPORT(loglik);
  REPORT(logdetpart);
  REPORT(randomContribution);
  REPORT(nu);

 
//  Rcout << "; logLik " <<  logLik << " nu " << nu << " loglik " << loglik;
  
  return -loglik;
}
