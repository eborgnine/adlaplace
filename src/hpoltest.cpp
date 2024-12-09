#include <TMB.hpp>

// Template function
template<class Type>
Type objective_function<Type>::operator() () {
  
  
  // Data inputs
  DATA_SPARSE_MATRIX(X);                   // Design matrix for fixed effects
  DATA_SPARSE_MATRIX(A);                   // Design matrix for random effects
  DATA_VECTOR(y);                   // Observed counts (Poisson responses)
  DATA_SPARSE_MATRIX(Q);            // Full precision matrix
  DATA_IVECTOR(gamma_split);         // Length of gamma for each city
  DATA_IMATRIX(cc_matrix);         // Length of gamma for each city
  int n_cc = cc_matrix.rows();
  int d_cc = cc_matrix.cols();
  
  
  // Parameter inputs
  PARAMETER_VECTOR(beta);                   // Fixed effects
  PARAMETER_VECTOR(gamma);                  // Concatenated random effects
  PARAMETER_VECTOR(theta);                  // Log-precision parameters (2J values for gamma0 and gamma)
  
  
  
  // Compute eta = log(lambda)
  Type nll = 0;
  vector<Type> eta_fixed = X * beta;
  vector<Type> eta_random = A * gamma;
  vector<Type> eta = eta_fixed + eta_random;
  
  // Compute negative log likelihood (multinomial -- case crossover)
  vector<Type> y_i(d_cc);
  vector<Type> eta_i(d_cc);
  Type lsa;
  vector<Type> prob_i(d_cc);
  for (int i = 0; i<n_cc; i++) {
    // reset
    y_i.setZero();
    eta_i.setZero();
    lsa = 0.0;
    prob_i.setZero();
    
    // y, eta and sum(eta) for ith colum of cc_matrix
    for(int j = 0; j<d_cc; j++) {
      if(cc_matrix(i,j) == 0) continue;
      y_i(j) = y(cc_matrix(i,j)-1);
      eta_i(j) = eta(cc_matrix(i,j)-1);
      lsa = logspace_add(lsa, eta(cc_matrix(i,j)-1));
    }
    
    // probabilities = exp(eta_i)/sum(exp(eta_i))
    for(int j = 0; j<d_cc; j++) {
      if(cc_matrix(i,j) == 0) continue;
      prob_i(j) = exp(logspace_sub(eta_i(j), lsa));
    }
    nll -= dmultinom(y_i, prob_i, 1);
  }
  // REPORT(nll);
  // Rcout << "ll : " << log_likelihood << "\n";





  // Random effects contribution
  vector<Type> exp_theta(A.cols());
  int idx_start = 0;
  int len_j = 0;
  for(int j = 0; j < gamma_split.size(); j++){
    len_j = gamma_split(j);
    // nll -= 0.5 * (len_j*theta(j) + log_dets(j)); // log-determinant part
    nll -= 0.5*len_j*theta(j); // log-determinant part
    for(int i = idx_start; i < idx_start + len_j; i++) exp_theta(i) = exp(theta(j)); // construct exp_theta for later
    idx_start += len_j;
  }

  // do sqrt(exp_theta) * gamma instead of exp_theta * Q
  gamma = gamma * sqrt(exp_theta);
  nll += 0.5*(gamma * (Q * gamma).col(0)).sum();
  
  return nll;
}
