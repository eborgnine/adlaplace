#include <TMB.hpp>
#define EIGEN_DONT_PARALLELIZE

//#define EVALCONSTANTS
#define LOGTWOPI 1.8378770664093454835606594728112352797227949472755668

template<class Type>
Type objective_function<Type>::operator() () {
  // Data inputs
  DATA_SPARSE_MATRIX(X);                   // Design matrix for fixed effects
  DATA_SPARSE_MATRIX(A);                   // Design matrix for random effects
  DATA_VECTOR(y);                          // Observed counts
  DATA_SPARSE_MATRIX(Q);                   // Precision matrix
  DATA_IVECTOR(gamma_split);               // Random effects grouping
  DATA_VECTOR(psd_scale);                  // Scaling factors
  DATA_IMATRIX(cc_matrix);                 // Case-control matrix
  
  // Parameters
  PARAMETER_VECTOR(beta);                  // Fixed effects
  PARAMETER_VECTOR(gamma);                 // Random effects
  PARAMETER_VECTOR(theta);                 // Variance parameters
  
  // Linear predictor (keep serial for AD)
  vector<Type> eta = X * beta + A * gamma;
  
  // Initialize parallel accumulator (thread-safe for AD)
  parallel_accumulator<Type> loglik_par(this);
  Type nu = theta(gamma_split.size());
  int n_cc = cc_matrix.rows();
  int d_cc = cc_matrix.cols();
  
  // Parallel region for case-control contributions
  PARALLEL_REGION {
    if (nu > 0) {
      // DIRICHLET-MULTINOMIAL parallelized
      Type logSqrtNu = log(nu)/2;
      Type oneOverSqrtNu = exp(-logSqrtNu);
      Type lgammaOneOverSqrtNu = lgamma(oneOverSqrtNu);
      
      for (int i = 0; i < n_cc; i++) {
        Type logSumMu = Type(-INFINITY);
        Type sumY = 0;
        
        // First pass (serial within parallel region)
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i,j);
          if (idx != 0) {
            idx -= 1;
            logSumMu = logspace_add(logSumMu, eta(idx));
            sumY += y(idx);
          }
        }
        
        Type contrib = lgammaOneOverSqrtNu - lgamma(oneOverSqrtNu + sumY);
#ifdef EVALCONSTANTS
        contrib += lgamma(1 + sumY);
#endif
        
        // Second pass
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i,j);
          if (idx != 0) {
            idx -= 1;
            Type muBarDivSqrtNu = exp(eta(idx) - logSumMu - logSqrtNu);
            Type yHere = y(idx);
            contrib += lgamma(yHere + muBarDivSqrtNu) - lgamma(muBarDivSqrtNu);
#ifdef EVALCONSTANTS
            contrib -= lgamma(yHere + 1);
#endif
          }
        }
        loglik_par += contrib; // Thread-safe accumulation
      }
    } else {
      // MULTINOMIAL parallelized
      for (int i = 0; i < n_cc; i++) {
        Type lsa = Type(-INFINITY);
#ifdef EVALCONSTANTS
        Type sumY = 0;
        Type sumGammaY_local = 0;
#endif
        
        // First pass
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i,j);
          if (idx != 0) {
            idx -= 1;
            lsa = logspace_add(lsa, eta(idx));
#ifdef EVALCONSTANTS
            sumY += y(idx);
            sumGammaY_local += lgamma(y(idx) + 1);
#endif
          }
        }
        
        // Second pass
        Type contrib = 0;
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i,j);
          if (idx != 0) {
            idx -= 1;
            contrib += y(idx) * (eta(idx) - lsa);
          }
        }
        
#ifdef EVALCONSTANTS
        contrib += lgamma(sumY + 1) - sumGammaY_local;
#endif
        loglik_par += contrib;
      }
    }
  } // End PARALLEL_REGION
  
  // Serial section for random effects (critical for AD)
  Type loglik = loglik_par; // Convert accumulator to regular Type
  
  // Random effects calculations (keep serial)
  vector<Type> sdTheta(A.cols());
  int idx_start = 0;
  for(int j = 0; j < gamma_split.size(); j++) {
    int len_j = gamma_split(j);
    for(int i = idx_start; i < idx_start + len_j; i++) {
      sdTheta(i) = theta(j)/psd_scale(j);
      loglik += log(sdTheta(i)); // Log-det part
    }
    idx_start += len_j;
  }
  
  gamma = gamma * sdTheta;
  Type randomContribution = -0.5*(gamma * (Q * gamma).col(0)).sum();
#ifdef EVALCONSTANTS
  randomContribution += (gamma.size()/2)*LOGTWOPI;
#endif
  
  loglik += randomContribution;
  ADREPORT(loglik);
  
  return -loglik;
}