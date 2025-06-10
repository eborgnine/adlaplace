#include "lgamma.hpp"
#include <Rcpp.h>

// [[Rcpp::export]]
double lgamma_forward_deriv(double x, int order) {
    if(order < 0 || order > 4)
        Rcpp::stop("Order must be between 0 and 4");

    CppAD::vector<double> tx(5, 0.0), ty(5, 0.0);
    tx[0] = x;
    tx[1] = 1.0; // For univariate, always seed tx[1]=1

    lgamma_ad_atomic.forward(0, CppAD::vector<bool>(1,true), 0, order, tx, ty);

    return ty[order];
}

// [[Rcpp::export]]
double lgamma_reverse_deriv(double x, int order) {
    if(order < 1 || order > 4)
        Rcpp::stop("Order must be between 1 and 4");

    CppAD::vector<double> tx(5, 0.0), ty(5, 0.0), px(5, 0.0), py(5, 0.0);
    tx[0] = x;
    tx[1] = 1.0;
    // First, run forward up to order-1 to fill Taylor coefficients
    lgamma_ad_atomic.forward(0, CppAD::vector<bool>(1,true), 0, order-1, tx, ty);

    py[order-1] = 1.0; // seed the adjoint for the previous derivative

    lgamma_ad_atomic.reverse(0, CppAD::vector<bool>(1,true), order-1, tx, ty, px, py);

    // px[0] is the order-th derivative
    return px[0];
}
