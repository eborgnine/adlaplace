#include "hpol.hpp"

#include <Rcpp.h>



// [[Rcpp::export]]
Rcpp::NumericMatrix logspaceadd_forward_deriv(Rcpp::NumericVector x, int order) {
    if(order < 0 || order > 4)
        Rcpp::stop("Order must be between 0 and 4");


    size_t n = x.size();
    Rcpp::NumericMatrix result(n, n);

    size_t q=order + 1;

    CppAD::vector<double> tx(q*n, 0.0), ty(q, 0.0);

    for (size_t i = 0; i < n; ++i)
        tx[i*q] = x[i];
    if(order == 1) {
        Rcpp::Rcout << "1";
        for (size_t i = 0; i < n; ++i) {
            tx[i*q+1] = 1;
            logspace_add_atomic.forward(0, CppAD::vector<bool>(1,true), 0, order, tx, ty);
            result(i,0) = ty[1];
            tx[i*q+1] = 0;
        }

    }
    if(order == 2) {
        Rcpp::Rcout << "\ndf2\n";
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
               tx[i*q+1] = 1;
               tx[j*q+1] = 1;
               logspace_add_atomic.forward(0, CppAD::vector<bool>(1,true), 0, order, tx, ty);
               result(i,j) = ty[2];
               tx[i*q+1] = 0;
               tx[j*q+1] = 0;
            }
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                result(i,j) -= result(i,i) + result(j,j);
                result(i,j) /= 2.0;
                result(j,i)  = result(i,j); 
            }
        }

    }

    return result;
}


// [[Rcpp::export]]
Rcpp::NumericVector logspaceadd_inbuilt_deriv(Rcpp::NumericVector x, int order) {
    size_t n = x.size();

    // Record with CppAD
    CppAD::vector<CppAD::AD<double>> ax(n);
    for (size_t i = 0; i < n; ++i)
        ax[i] = x[i];
    CppAD::Independent(ax);

    CppAD::vector<CppAD::AD<double>> ay(1);
    ay[0] = logspace_add_n(ax);

    CppAD::ADFun<double> f(ax, ay);

    // Compute gradient
    std::vector<double> x0(n);
    for (size_t i = 0; i < n; ++i) x0[i] = x[i];
  if(order == 0) {
        double val = f.Forward(0, x0)[0];
        return Rcpp::NumericVector::create(val);
    }
    if(order == 1) {
        std::vector<double> grad = f.Jacobian(x0);
        return Rcpp::wrap(grad);
    }

 if(order == 2) {
        // Compute Hessian of f(x) wrt x
        std::vector<double> w(1, 1.0); // seed for output
        std::vector<double> hess = f.Hessian(x0, 0);
        // Return as matrix
        Rcpp::NumericMatrix mat(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                mat(i, j) = hess[i * n + j];
            }
        }
        return mat;
    }
 if(order == -2) {
        // Compute Hessian of f(x) wrt x
    std::vector<double> u1(n, 0), u2(n,0);
    Rcpp::NumericMatrix mat(n, n);
    for (size_t i = 0; i < n; ++i){
        u1[i]=1;
        f.Forward(1, u1);
        auto y2 = f.Forward(2, u2);
        mat(i,i) = y2[0];
        u1[i]=0;
        for (size_t j = 0; j < i; ++j){
           u1[i]=1;
           u1[j]=1;
           f.Forward(1, u1);
           y2 = f.Forward(2, u2);
           mat(i,j) = (y2[0]- mat(i,i) - mat(j,j))/2;
           mat(j,i) = mat(i,j);
           u1[i]=0;
           u1[j]=0;
        }
    }
    return mat;
} 

if(order == 3) {
    Rcpp::NumericVector third(n * n * n);
    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < n; ++j)
            for(size_t k = 0; k < n; ++k) {
                std::vector<double> dx(n, 0.0), ddx(n, 0.0), dddx(n, 0.0);
                dx[i] = 1.0; ddx[j] = 1.0; dddx[k] = 1.0;
                // Prepare tx for Forward(3, ...)
                std::vector<double> tx(n*4, 0.0);
                for(size_t l = 0; l < n; ++l) {
                    tx[l + 0*n] = x0[l];
                    tx[l + 1*n] = dx[l];
                    tx[l + 2*n] = ddx[l];
                    tx[l + 3*n] = dddx[l];
                }
                auto y = f.Forward(3, tx);
                third[i + n*j + n*n*k] = y[0] * 6.0; // third derivative
            }
    Rcpp::NumericVector thirddim(3);
    thirddim[0] = n;
    thirddim[1] = n;
    thirddim[2] = n;

    third.attr("dim") = thirddim;
    return third;
    }


    return Rcpp::wrap(0);
}



CppAD::vector< CppAD::AD<double> > testAddGamma(  
    CppAD::vector<CppAD::AD<double>> ad_params) {


    ad_params[0] = ad_params[0]*ad_params[0];
    ad_params[1] = ad_params[1]*ad_params[0];
    ad_params[2] = lgamma_ad(ad_params[2]);

    CppAD::vector< CppAD::AD<double> > result(1);
    result[0] = logspace_add_n(ad_params);
    return result;
}


// [[Rcpp::export]]
Rcpp::List testAddGammaR(
  Rcpp::NumericVector parameters
  ) {

    size_t Nparams = parameters.size();
    CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
    for (size_t D = 0; D < Nparams; D++) {
        ad_params[D] = parameters[D];  // Initialize CppAD variables
    }
    CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

    auto y = testAddGamma(ad_params);
    CppAD::ADFun<double> f(ad_params, y);

    std::vector<double> x_val(Nparams);
    for (size_t i = 0; i < Nparams; ++i) x_val[i] = parameters[i];

    std::vector<double> y_val(1);
    y_val = f.Forward(0, x_val);
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("value") = y_val[0]
    );

  std::vector<double> grad = f.Jacobian(x_val);

  result["grad"] = grad;

  std::vector<double> hess = f.Hessian(x_val, 0);
        // Return as matrix
        Rcpp::NumericMatrix mat(Nparams, Nparams);
        for (size_t i = 0; i < Nparams; ++i) {
            for (size_t j = 0; j < Nparams; ++j) {
                mat(i, j) = hess[i * Nparams + j];
            }
        }
  result["hessian"] = mat;

  return result;
}
