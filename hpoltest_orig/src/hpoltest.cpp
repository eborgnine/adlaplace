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
  DATA_VECTOR(psd_scale);                  // Scaling factors for predictive standard deviations
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
  int n_cc = cc_matrix.rows();            // Number of case-control groups
  int d_cc = cc_matrix.cols();             // Max group size
  
  PARALLEL_REGION {
    if (dirichlet) {
      // Precompute Dirichlet constants (thread-safe)
      Type nu = theta.tail(1)(0); 
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
        
        Type contrib = lgammaOneOverSqrtNu - lgamma(oneOverSqrtNu + sumY);
#ifdef EVALCONSTANTS
        contrib += lgamma(sumY + 1);
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
  
  Type loglik = loglik_par;  // Convert parallel accumulator to serial

  
  // ======================
  // 5. Random Effects Penalty (Serial)
  // ======================
  
  // SD with PSD
  
  vector<Type> thetaSub = theta.head(gamma_split.size());
  if (psd_scale.size() > 0) {
    vector<Type> psd_scale_cast = psd_scale.head(thetaSub.size()).template cast<Type>();
//    thetaSub = thetaSub / psd_scale_cast;  // Apply scaling
  }  
  
  Type randomContribution = -(gamma_split.template cast<Type>() * log(thetaSub)).sum();
  
  
  // Scale random effects by theta
  vector<Type> scaled_gamma(gamma.size());
  int idx = 0;
  for (int j = 0; j < gamma_split.size(); j++) {
    for (int k = 0; k < gamma_split(j); k++) {
      scaled_gamma(idx) = gamma(idx) / thetaSub(j);  // Scale by group theta
      idx++;
    }
  }
  randomContribution -= 0.5 * (scaled_gamma * (Q * scaled_gamma)).sum();

#ifdef EVALCONSTANTS
  randomContribution -= 0.5*gamma.size() * LOGTWOPI;
  randomContribution += 0.5*logdet(Q);
#endif
  
  loglik += randomContribution;

  REPORT(gamma);
  REPORT(eta);
  REPORT(beta);
  REPORT(theta);
  REPORT(thetaSub);
  REPORT(randomContribution);
  REPORT(loglik);

  return -loglik;  // Minimize negative log-likelihood
}