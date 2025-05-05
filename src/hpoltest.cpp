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
  // final parameter is nu
  Type nu = theta(gamma_split.size()+1);

  // Compute eta = log(lambda)
  Type nll = 0;
  vector<Type> eta(y.size());
  eta.setZero();
  if(X.cols() > 0) eta += X * beta;
  if(A.cols() > 0) eta += A * gamma;

  // Compute negative log likelihood (multinomial -- case crossover)
  Type lsa;

  if(nu > 0) { // dirichelet multinom

    int logSqrtNu = log(nu)/2;
    int oneOverSqrtNu = exp(-logOneOverSqrtNu);
    int gammaNu = lgamma(oneOverSqrtNu);
    int sumY, gammaNuSumY, muBar, yHere;

    for (int i = 0; i<n_cc; i++) {
      lsa = Type(-INFINITY);
      sumY = 0L;
      for(int j = 0; j<d_cc; j++)
        lsa = Type(-INFINITY);
        if(cc_matrix(i,j) != 0) {
          lsa = logspace_add(lsa, eta(cc_matrix(i,j)-1));
          sumY += y(cc_matrix(i,j)-1) 
        }
        gammaNuSumY = gammaNu + lgamma(1+sumY) - lgamma(oneOverSqrtNu+sumY);
        lsa = lsa + logSqrtNu;

    // etaBar = (1/sqrt(nu)) * exp(eta_ij)/sum(exp(eta_i)) = exp(eta_ij - lsa)
        for(int j = 0; j<d_cc; j++) 
          if(cc_matrix(i,j) != 0) {
            muBar = exp(eta(cc_matrix(i,j)-1) - lsa);
            yHere = y(cc_matrix(i,j)-1);
            nll -= lgamma(yHere + muBar) - lgamma(yHere + 1) - lgamma(muBar);
          }
        }
  } else { // multinom
    for (int i = 0; i<n_cc; i++) {
      lsa = Type(-INFINITY);
      for(int j = 0; j<d_cc; j++) 
        if(cc_matrix(i,j) != 0) lsa = logspace_add(lsa, eta(cc_matrix(i,j)-1));

    // - y_ij * log(p_ij), where p_ij = exp(eta_ij)/sum(exp(eta_i)) = exp(eta_ij - lsa)
      for(int j = 0; j<d_cc; j++) 
        if(cc_matrix(i,j) != 0) nll -= y(cc_matrix(i,j)-1) * (eta(cc_matrix(i,j)-1) - lsa);
    }
  } // end multinom
  // REPORT(nll);





  // Random effects contribution
  vector<Type> exp_theta(A.cols());
  int idx_start = 0;
  int len_j = 0;
  for(int j = 0; j < gamma_split.size(); j++){
    len_j = gamma_split(j);
    // nll -= 0.5 * (len_j*theta(j) + log_dets(j)); // log-determinant part
    nll -= 0.5*len_j*theta(j);
    for(int i = idx_start; i < idx_start + len_j; i++) exp_theta(i) = exp(theta(j)); // construct exp_theta for later
      idx_start += len_j;
  }

  // do sqrt(exp_theta) * gamma instead of exp_theta * Q
  gamma = gamma * sqrt(exp_theta);
  nll += 0.5*(gamma * (Q * gamma).col(0)).sum();

  return nll;
}
