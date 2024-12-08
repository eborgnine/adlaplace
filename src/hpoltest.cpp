#include <TMB.hpp>

// Template function
template<class Type>
Type objective_function<Type>::operator() () {
  
  
  // Data inputs
  DATA_SPARSE_MATRIX(X);                   // Design matrix for fixed effects
  DATA_SPARSE_MATRIX(A);                   // Design matrix for random effects
  DATA_VECTOR(y);                   // Observed counts (Poisson responses)
  DATA_SPARSE_MATRIX(Q);            // Full precision matrix
  DATA_VECTOR(log_dets);              // log determinant of exp(theta) * Q
  DATA_IVECTOR(gamma_split);         // Length of gamma for each city
  DATA_IVECTOR(theta_id);            // Length of gamma for each city
  DATA_INTEGER(nc)              // Number of cities
  DATA_INTEGER(np)              // Number of exposures
    
  
  // Parameter inputs
  PARAMETER_VECTOR(beta);                   // Fixed effects
  PARAMETER_VECTOR(gamma);                  // Concatenated random effects
  PARAMETER_VECTOR(theta);                  // Log-precision parameters (2J values for gamma0 and gamma)
  
  
  
  // Initialize negative log-likelihood
  Type nll = 0;
  vector<Type> eta_fixed = X * beta;
  vector<Type> eta_random = A * gamma;
  vector<Type> log_lambda = eta_fixed + eta_random;
  for (int i = 0; i < y.size(); i++) {
    nll -= dpois(y(i), exp(log_lambda(i)), true);
  }



  // Random effects contribution
  vector<Type> exp_theta(A.cols());
  int idx_start = 0;
  int len_j = 0;
  int log_det_id = 0;
  for(int j = 0; j < gamma_split.size(); j++){
    len_j = gamma_split(j);
    for(int i = idx_start; i < idx_start + gamma_split(j); i++){
      exp_theta(i) = exp(theta(theta_id(j)));                     // construct exp_theta for later
    }
    nll -= 0.5 * (gamma_split(j)*theta(theta_id(j)) + log_dets(log_det_id)); // log-determinant part
    log_det_id += 1;
    if(log_det_id == np) log_det_id = 0;
    idx_start += len_j;
  }

  // do sqrt(exp_theta) * gamma instead of exp_theta * Q
  gamma = gamma * sqrt(exp_theta);
  int d = Q.cols();
  vector<Type> gamma_i(d);
  for (int i = 0; i < (nc+1); i++) {
    gamma_i = gamma.segment(i*d, d);
    nll += 0.5 * (gamma_i * (Q * gamma_i).col(0)).sum();
  }

  return nll;
}
