#include <TMB.hpp>
#define EIGEN_DONT_PARALLELIZE

// Constants (consider moving to R-side if they change often)
#define LOGTWOPI 1.8378770664093454835606594728112352797227949472755668

template<class Type>
Type objective_function<Type>::operator() () {
  // ======================
  // 1. Data Inputs
  // ======================
  DATA_SPARSE_MATRIX(X);                   // Fixed effects design matrix
  DATA_SPARSE_MATRIX(A);                   // Random effects design matrix
  DATA_VECTOR(y);                          // Observed counts
  DATA_SPARSE_MATRIX(Q);                   // Precision matrix
  DATA_IVECTOR(gamma_split);               // Random effects grouping
  DATA_VECTOR(psd_scale);                  // Scaling factors (unused in current code)
  DATA_IMATRIX(cc_matrix);                 // Case-control matrix
  DATA_INTEGER(dirichlet);                 // 0=Multinomial, 1=Dirichlet-Multinomial
  
  // ======================
  // 2. Parameters
  // ======================
  PARAMETER_VECTOR(beta);                  // Fixed effects
  PARAMETER_VECTOR(gamma);                 // Random effects
  PARAMETER_VECTOR(theta);                 // Variance parameters
  
  // ======================
  // 3. Linear Predictor
  // ======================
  vector<Type> eta = X * beta + A * gamma;
  
  if (theta(theta.size()-1) < 0) {
    Rf_error("theta must be positive!", theta(theta.size()-1));
  }
  
  // ======================
  // 4. Parallel Likelihood Calculation
  // ======================
  parallel_accumulator<Type> loglik_par(this);
  Type nu = theta(gamma_split.size());     // Dirichlet dispersion parameter
  int n_cc = cc_matrix.rows();            // Number of case-control groups
  int d_cc = cc_matrix.cols();             // Max group size
  
  PARALLEL_REGION {
    if (dirichlet) {
      // Precompute Dirichlet constants (thread-safe)
      Type logSqrtNu = log(nu) / 2;
      Type oneOverSqrtNu = exp(-logSqrtNu);
      Type lgammaOneOverSqrtNu = lgamma(oneOverSqrtNu);
      
      for (int i = 0; i < n_cc; i++) {
        Type logSumMu = Type(-INFINITY);
        Type sumY = 0;
        
        // First pass: Compute logSumMu and sumY
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i, j);
          if (idx != 0) {
            idx -= 1;  // Convert to 0-based index
            logSumMu = logspace_add(logSumMu, eta(idx));
            sumY += y(idx);
          }
        }
        
        // Dirichlet-Multinomial contribution
        Type contrib = lgammaOneOverSqrtNu - lgamma(oneOverSqrtNu + sumY);
#ifdef EVALCONSTANTS
        contrib += lgamma(1 + sumY);
#endif
        
        // Second pass: Add per-observation terms
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i, j);
          if (idx != 0) {
            idx -= 1;
            Type muBarDivSqrtNu = exp(eta(idx) - logSumMu - logSqrtNu);
            contrib += lgamma(y(idx) + muBarDivSqrtNu) - lgamma(muBarDivSqrtNu);
#ifdef EVALCONSTANTS
            contrib -= lgamma(y(idx) + 1);
#endif
          }
        }
        loglik_par += contrib;
      }
    } else {
      // Multinomial likelihood (parallelized)
      for (int i = 0; i < n_cc; i++) {
        Type logSumExpEta = Type(-INFINITY);
#ifdef EVALCONSTANTS
        Type sumY = 0;
        Type sumLgammaY = 0;
#endif
        
        // First pass: Compute logSumExpEta
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i, j);
          if (idx != 0) {
            idx -= 1;
            logSumExpEta = logspace_add(logSumExpEta, eta(idx));
#ifdef EVALCONSTANTS
            sumY += y(idx);
            sumLgammaY += lgamma(y(idx) + 1);
#endif
          }
        }
        
        // Second pass: Compute multinomial log-likelihood
        Type contrib = 0;
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i, j);
          if (idx != 0) {
            idx -= 1;
            contrib += y(idx) * (eta(idx) - logSumExpEta);
          }
        }
#ifdef EVALCONSTANTS
        contrib += lgamma(sumY + 1) - sumLgammaY;
#endif
        loglik_par += contrib;
      }
    }
  }  // End PARALLEL_REGION
  
  // ======================
  // 5. Random Effects Penalty (Serial)
  // ======================
  Type loglik = loglik_par;  // Convert parallel accumulator to serial
  
  // Scale random effects by theta
  vector<Type> sdTheta(A.cols());
  int idx_start = 0;
  for (int j = 0; j < gamma_split.size(); j++) {
    int len_j = gamma_split(j);
    Type theta_j = theta(j);  // Reuse for all elements in group j
    for (int i = idx_start; i < idx_start + len_j; i++) {
      sdTheta(i) = theta_j;  // Note: psd_scale unused (commented out)
      loglik += log(sdTheta(i));  // Log-det of transformation
    }
    idx_start += len_j;
  }
  
  gamma = gamma * sdTheta;
  Type randomContribution = -0.5 * (gamma * (Q * gamma).col(0)).sum();
#ifdef EVALCONSTANTS
  randomContribution += (gamma.size() / 2) * LOGTWOPI;
#endif
  
  loglik += randomContribution;
  REPORT(eta);
  REPORT(beta);
  REPORT(theta);
  REPORT(randomContribution);
  REPORT(loglik);

  return -loglik;  // Minimize negative log-likelihood
}