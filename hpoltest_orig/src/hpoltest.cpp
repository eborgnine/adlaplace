

#include <cmath> 
#include <Rcpp.h>
#include <TMB.hpp>
#include "lgamma_ad.hpp"

#define EIGEN_DONT_PARALLELIZE

// Constants (consider moving to R-side if they change often)
#define LOGTWOPI 1.8378770664093454835606594728112352797227949472755668



/*
 * LGAMMA
 */

// Define the global atomic function instance
atomic_lgamma_ad lgamma_ad_atomic("lgamma_ad_atomic");

// Constructor implementation
atomic_lgamma_ad::atomic_lgamma_ad(const std::string& name)
    : CppAD::atomic_base<double>(name) {}



// Forward mode
bool atomic_lgamma_ad::forward(
    size_t p,
    size_t q,
    const CppAD::vector<bool>& vx,
    CppAD::vector<bool>& vy,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
) {

     double deriv[4], tx1mult=0.0;

    // Always fill derivatives up to max requested order!
    deriv[0] = R::psigamma(tx[0], 0); // digamma
    if (q >= 1)
        deriv[1] = R::psigamma(tx[0], 1); // trigamma
    if (q >= 2)
        deriv[2] = R::psigamma(tx[0], 2);
    if (q >= 3)
        deriv[3] = R::psigamma(tx[0], 3);

    // Zero order (function value)
    if (p <= 0)
        ty[0] = std::lgamma(tx[0]);

    if (q >= 1 && p <= 1)
        ty[1] = deriv[0] * tx[1];

    if (q >= 2 && p <= 2) {
        tx1mult = tx[1] * tx[1];
        ty[2] = deriv[0] * tx[2] + deriv[1] * tx1mult;
    }

    if (q >= 3 && p <= 3) {
        tx1mult *= tx[1];
        ty[3] = deriv[0] * tx[3]
              + 3 * deriv[1] * tx[1] * tx[2]
              + deriv[2] * tx1mult;
    }

    if (q >= 4 && p <= 4) {
        tx1mult *= tx[1];
        ty[4] = deriv[0] * tx[4]
              + 4 * deriv[1] * tx[1] * tx[3]
              + 3 * deriv[1] * tx[2] * tx[2]
              + 6 * deriv[2] * tx[1] * tx[1] * tx[2]
              + deriv[3] * tx1mult;
    }

    return true;
}

// Reverse mode
bool atomic_lgamma_ad::reverse(
    size_t q,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
) {

    for(size_t i = 0; i < px.size(); ++i)
      px[i] = 0.0;

    double x0 = tx[0];
    double deriv[4];
    double tx1sq, tx1cubed;


    // Order 0
    if (q >= 0) {
      deriv[0] = R::psigamma(x0, 0);
      px[0] = py[0] * deriv[0];
    }

    // Order 1
    if (q >= 1) {
      deriv[1] = R::psigamma(x0, 1);
      px[0] += py[1] * deriv[1] * tx[1];
      px[1] = py[1] * deriv[0];
    }

    // Order 2
    if (q >= 2) {
      deriv[2] = R::psigamma(x0, 2);      
      tx1sq = tx[1]*tx[1];
      px[0] += py[2] * (deriv[2] * tx1sq + deriv[1] * tx[2]);
      px[1] += py[2] * 2 * deriv[1] * tx[1];
      px[2] = py[2] * deriv[0];
    }

  // Order 3
    if (q >= 3) {
        tx1cubed = tx[1]*tx1sq;
        deriv[3] = R::psigamma(x0, 3);      
        px[0] += py[3] * 
          (deriv[1] * tx[3] + 
            3 * deriv[2] * tx[1] * tx[2] + 
            deriv[3] * tx1cubed);
        px[1] += py[3] * (
            3 * deriv[1] * tx[2] + 
            3 * deriv[2] * tx1sq);
        px[2] += py[3] * (3 * deriv[0] * tx[1]);
        px[3] += py[3] * deriv[0];
    }
    return true;
  }




/* 
* logspace add
*/ 

template<class Type>
Type logspace_add_n(const CppAD::vector<Type>& x) {

    auto n = x.size();
    Type max_x = x[0];
    for (size_t i = 1; i < n; ++i) {
        if (x[i] > max_x) {
            max_x = x[i];
        }
    }

    Type sum_exp = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum_exp += exp(x[i] - max_x);
    }
    return max_x + log(sum_exp);
}

/*
* MODEL
*/



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

#ifdef DEBUG
  Rcpp::Rcout << "sum log etas\n";
#endif


  if (theta(theta.size()-1) < 0) {
    Rcpp::Rcout << "theta must be positive! " <<
     theta(theta.size()-1) << "\n";
  }
  
  // ======================
  // 4. Parallel Likelihood Calculation
  // ======================
  parallel_accumulator<Type> loglik_par(this);
  int n_cc = cc_matrix.rows();            // Number of case-control groups
  int d_cc = cc_matrix.cols();             // Max group size
  
  PARALLEL_REGION {
      
    for (int i = 0; i < n_cc; i++) {
        Type sumY = 0;
        int NinStrata = 0; // count number in this strata
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i, j);
          if (idx != 0) {
            NinStrata ++;
          }
        }
        CppAD::vector<Type> etaHere(NinStrata);
        for (int j = 0; j < d_cc; j++) {
          int idx = cc_matrix(i, j);
          if (idx != 0) {
            etaHere[j] = eta[idx];
            sumY += y(idx);
          }
        }
        Type logSumMu = logspace_add_n(etaHere);
        Type contrib=0;    

        if (dirichlet) {
          Type nu = theta.tail(1)(0); 
          Type logSqrtNu = log(nu) / 2;
          Type oneOverSqrtNu = exp(-logSqrtNu);
          Type lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu);

          contrib += lgammaOneOverSqrtNu - lgamma_ad(oneOverSqrtNu + sumY);
          for (int j = 0; j < d_cc; j++) {
            int idx = cc_matrix(i, j);
            if (idx != 0) {
              idx -= 1;
              Type muBarDivSqrtNu = exp(eta[idx] - logSumMu - logSqrtNu);
              contrib += lgamma_ad(y(idx) + muBarDivSqrtNu) - lgamma_ad(muBarDivSqrtNu);
            }
          }
          loglik_par += contrib;
        } else {
      // Multinomial likelihood (parallelized)
          for (int i = 0; i < n_cc; i++) {        
            for (int j = 0; j < d_cc; j++) {
              int idx = cc_matrix(i, j);
              if (idx != 0) {
                idx -= 1;
                contrib += y(idx) * (eta(idx) - logSumMu);
              }
            }
            loglik_par += contrib;
          }
        } 
    } // for strata
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

  
  loglik += randomContribution;

  return -loglik;  // Minimize negative log-likelihood
}