#include "logspaceadd.hpp"
#include <Rcpp.h>


// Constructor implementation
atomic_logspace_add::atomic_logspace_add(const std::string& name, size_t n_)
    : CppAD::atomic_four<double>(name), n(n_) {}


bool atomic_logspace_add::for_type(
    size_t call_id,
    const CppAD::vector<CppAD::ad_type_enum>& type_x,
    CppAD::vector<CppAD::ad_type_enum>& type_y
  ) {
    // For log-sum-exp with n inputs, output is scalar
    assert(type_y.size() == 1);
    // Take max type of any input as output type
    type_y[0] = type_x[0];
    for(size_t i = 1; i < type_x.size(); ++i)
        type_y[0] = std::max(type_y[0], type_x[i]);
    return true;
}

bool atomic_logspace_add::forward(
    size_t call_id,
    const CppAD::vector<bool>& select_y,
    size_t order_low,
    size_t order_up,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
  )  {
//    size_t n_order = order_up + 1;
    size_t n = this->n; // number of inputs
    size_t q = order_up + 1;


    std::vector<double> x0(n), dx(n), dx2(n), ddx(n);
    for(size_t i=0; i<n; ++i) {
        x0[i]  = tx[i * q + 0];
    }
    auto max_iter = std::max_element(x0.begin(), x0.end());
    size_t max_idx = std::distance(x0.begin(), max_iter);
    double max_x = *max_iter;

  
    // Numerically stable log-sum-exp
    double sum_exp = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (x0[i] == max_x) continue; // skip max for stability
        sum_exp += std::exp(x0[i] - max_x);
    }
    double logsumexp = max_x + std::log1p(sum_exp);


    if (order_low <= 0)
        ty[0] = logsumexp;


    // Softmax weights (include max index)
    std::vector<double> w(n);

 // First order
    if (order_up >= 1) {
        double t1 = 0.0;
        for(size_t i = 0; i < n; ++i){
            dx[i] = tx[i * q + 1];
            w[i] = std::exp(x0[i] - logsumexp);
            t1 += w[i] * dx[i];
        }
        ty[1] = t1;
    }

    // Second order (Hessian applied to direction)
    if (order_up >= 2) {

        ty[2] =0;

        for (size_t j = 0; j < n; ++j) {
            dx2[j] = tx[j * q + 2]; // second direction
            ty[2] += exp(x0[j]) /  exp(logsumexp) * dx2[j];
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ty[2] +=  ( // hess
                    (i == j) * exp(x0[i])/exp(logsumexp) 
                    - exp(x0[i] + x0[j]) / exp(2*logsumexp)
                ) * dx[i] *dx[j];
            }
        }
    } // order 2
   
    if (order_up >= 3) {
          // We'll need y1 for the higher orders
    double y1 = 0.0;
    std::vector<double> wi(n);
    for (size_t i = 0; i < n; ++i) {
        wi[i] = std::exp(tx[i*q + 0] - max_x) / sum_exp;
        y1 += wi[i] * tx[i*q + 1];
    }
     // Compute y3 according to the formula:
    double y3 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double delta1 = tx[i*q + 1] - y1;
        y3 += wi[i] * tx[i*q + 3];
        y3 += 3.0 * wi[i] * delta1 * tx[i*q + 2];
        y3 += wi[i] * delta1 * delta1 * delta1;
    }

    ty[3] = y3;
    }

    return true;
  }

bool atomic_logspace_add::reverse(
    size_t call_id,
    const CppAD::vector<bool>& select_x,
    size_t order_up,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
)  {
    size_t q = order_up + 1;
    size_t n = this->n; // number of inputs

    assert(tx.size() >= 2 * q);
    assert(ty.size() >= 1 * q);
    assert(px.size() >= 2 * q);
    assert(py.size() >= 1 * q);

    std::fill(px.begin(), px.end(), 0.0);

    // Zeroth order (function value)
    std::vector<double> x0(n), x1(n), x2(n),xDivSum(n);
    for(size_t i = 0; i < n; ++i) {
        x0[i] = tx[i*q];  
        x1[i] = tx[i*q+1];  
        x2[i] = tx[i*q + 2];
    }

    auto max_iter = std::max_element(x0.begin(), x0.end());
    size_t max_idx = std::distance(x0.begin(), max_iter);
    double max_x = *max_iter;

    // Numerically stable log-sum-exp
    double sum_exp = 0.0;
    std::vector<double> expDiff(n);
    for (size_t i = 0; i < n; ++i) {
        expDiff[i] = std::exp(x0[i] - max_x);
        if (x0[i] == max_x) continue; // skip max for stability
        sum_exp += expDiff[i];
    }
    double logsumexp = max_x + std::log1p(sum_exp), expMaxX = exp(max_x);
    for (size_t i = 0; i < n; ++i) {
        xDivSum[i] = std::exp(x0[i]-logsumexp);
    }

// Compute S and S1
    double S = exp(logsumexp), S1 = 0.0, S2=0.0;
    for(size_t i = 0; i < n; ++i) {
        S1 += expDiff[i] * x1[i];
        S2 += expDiff[i] * x2[i];
    }
    double logS1 = max_x + log(S1), logS2 = max_x + log(S2);
    S1 = exp(logS1);
    S2 = exp(logS2);
    double logS1divS = logS1 - logsumexp;
    double S1divS = exp(logS1divS), S1divSsq = exp(2*logS1divS);

// order 0
    if(order_up >= 0) {

    for(size_t i = 0; i < n; ++i)
        px[i*q] = py[0] * xDivSum[i];
    }
// order 1
    if(order_up >= 1) {
    for(size_t i = 0; i < n; ++i) {
        px[q*i] += py[1] * xDivSum[i] * (x1[i] - S1divS);
        px[q * i + 1] = py[1] * xDivSum[i];
    }
}

    if(order_up >= 2) {
        for(size_t i = 0; i < n; ++i) {
            double vi = x2[i] - 2 * S1divS * x1[i] + 
                2 *  S1divSsq - S2/S;
            px[i*q + 2] = py[2] * xDivSum[i];
            px[i*q + 0] += py[2] * xDivSum[i] * vi;

            // The "mixed" second-order terms, which couple x1 for i,j
            for (size_t j = 0; j < n; ++j) {
                px[i*q + 0] += py[2] * xDivSum[i] *
                    (x1[i] - S1divS) * (i==j ? 1 : 0) * x1[j];
                px[i*q + 0] -= py[2] * xDivSum[i] * (x1[i] - S1divS) * xDivSum[j] * x1[j];
            }
        }
    }
    return true;
}


