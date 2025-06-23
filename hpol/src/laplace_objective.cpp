#include "hpol.hpp"
#include <RcppEigen.h>

// Taping the inner gradient (AD<double>)
CppAD::ADFun<double> tape_inner_gradient(
    const std::vector<double>& gamma_hat_numeric,
    const std::vector<double>& beta_theta_numeric,  
    const Rcpp::List& data,
    const Rcpp::List& config
) {
    size_t Ngamma = gamma_hat_numeric.size();
    size_t NbetaTheta = beta_theta_numeric.size();

    CppAD::vector<CppAD::AD<double>> gamma_hat_AD(Ngamma);
    for (size_t i = 0; i < Ngamma; ++i)
        gamma_hat_AD[i] = gamma_hat_numeric[i];

    CppAD::vector<CppAD::AD<double>> beta_theta_AD(NbetaTheta);
    for (size_t i = 0; i < NbetaTheta; ++i)
        beta_theta_AD[i] = beta_theta_numeric[i];

    CppAD::Independent(beta_theta_AD);

    auto grad_inner = [&](const CppAD::vector<CppAD::AD<double>>& gamma_inner) -> CppAD::vector<CppAD::AD<double>> {
        CppAD::vector<CppAD::AD<double>> full_vec(NbetaTheta + Ngamma);
        for (size_t i = 0; i < NbetaTheta; ++i)
            full_vec[i] = beta_theta_AD[i];
        for (size_t i = 0; i < Ngamma; ++i)
            full_vec[NbetaTheta + i] = gamma_inner[i];

        return objectiveFunctionGradientGamma(full_vec, data, config); 
    };

    // Gradient at gamma_hat (AD<double>)
    CppAD::vector<CppAD::AD<double>> grad_at_gamma_hat = grad_inner(gamma_hat_AD);

    CppAD::ADFun<double> inner_tape(beta_theta_AD, grad_at_gamma_hat);
    return inner_tape; 
}

// Outer tape: Derivative of gamma_hat w.r.t beta_theta (AD<AD<double>>)
CppAD::ADFun<CppAD::AD<double>> tape_outer_jacobian(
    CppAD::ADFun<double>& inner_tape,
    const std::vector<double>& beta_theta_numeric
) {
    size_t NbetaTheta = beta_theta_numeric.size();
    CppAD::vector<CppAD::AD<CppAD::AD<double>>> beta_theta_outer(NbetaTheta);
    for(size_t i = 0; i < NbetaTheta; ++i) 
        beta_theta_outer[i] = beta_theta_numeric[i];

    CppAD::Independent(beta_theta_outer);

    std::vector<CppAD::AD<CppAD::AD<double>>> beta_theta_inner(NbetaTheta);
    for (size_t i = 0; i < NbetaTheta; ++i)
        beta_theta_inner[i] = CppAD::Value(beta_theta_outer[i]);

    // Evaluate inner_tape at beta_theta_inner
    std::vector<CppAD::AD<CppAD::AD<double>>> inner_result = inner_tape.Forward(0, beta_theta_inner);

    CppAD::ADFun<CppAD::AD<double>> outer_tape(beta_theta_outer, inner_result);
    return outer_tape;
}

// R interface
// [[Rcpp::export]]
Rcpp::List loglik_laplace(
    Rcpp::NumericVector beta_theta,
    Rcpp::NumericVector gamma_hat_numeric,
    Rcpp::List data,
    Rcpp::List config
) {
    std::vector<double> beta_theta_numeric(beta_theta.begin(), beta_theta.end());
    std::vector<double> gamma_hat_std(gamma_hat_numeric.begin(), gamma_hat_numeric.end());

    // Inner tape: gradient w.r.t gamma
    auto inner_fun = tape_inner_gradient(gamma_hat_std, beta_theta_numeric, data, config);

    // Outer tape: jacobian of gamma_hat w.r.t beta_theta
    auto outer_fun = tape_outer_jacobian(inner_fun, beta_theta_numeric);

    // Numeric evaluation of Jacobian at beta_theta_numeric
    std::vector<double> jacobian_vec = outer_fun.Jacobian(beta_theta_numeric);

    size_t Ngamma = gamma_hat_numeric.size();
    size_t NbetaTheta = beta_theta_numeric.size();
    Rcpp::NumericMatrix jacobian(Ngamma, NbetaTheta);

    for (size_t row = 0; row < Ngamma; ++row) {
        for (size_t col = 0; col < NbetaTheta; ++col) {
            jacobian(row, col) = jacobian_vec[row * NbetaTheta + col];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("jacobian") = jacobian
    );
}
