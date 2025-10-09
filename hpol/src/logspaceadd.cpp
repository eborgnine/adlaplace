#ifdef UNDEF

#include "logspaceadd.hpp"
#include <Rcpp.h>




// Constructor implementation
atomic_logspace_add::atomic_logspace_add(const std::string& name)
    : CppAD::atomic_four<double>(name) {}


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
    size_t q = order_up + 1;
    size_t n = tx.size()/q;// number of inputs


    std::vector<double> x0(n), dx(n), dx2(n), ddx(n);
    for(size_t i=0; i<n; ++i) {
        x0[i]  = tx[i * q + 0];
    }
    auto max_iter = std::max_element(x0.begin(), x0.end());
    double max_x = *max_iter;
    size_t max_idx = std::distance(x0.begin(), max_iter);

//    Rcpp::Rcout << "max_x " << max_x << "\n";
    // Numerically stable log-sum-exp
    double sum_exp = 0.0;
    for (size_t i = 0; i < n; ++i) {
        if (i == max_idx) continue; // skip max for stability
        sum_exp += std::exp(x0[i] - max_x);
    }
    double logsumexp = max_x + std::log1p(sum_exp);
//    Rcpp::Rcout << "sum_exp " << sum_exp << " logsumexp " << logsumexp << "\n";


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
            ty[2] += exp(x0[j] - logsumexp) * dx2[j];
        }
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                ty[2] +=  ( // hess
                    (i == j) * exp(x0[i]-logsumexp) 
                    - exp(x0[i] + x0[j] - 2*logsumexp)
                ) * dx[i] *dx[j];
            }
        }
    } // order 2
   
    if (order_up >= 3) {
              Rcpp::Rcout << "logspaceadd forward 3 havent tested this\n";
          // We'll need y1 for the higher orders
    double y1 = 0.0;
    std::vector<double> wi(n);
    for (size_t i = 0; i < n; ++i) {
        wi[i] = std::exp(tx[i*q + 0] - max_x) / sum_exp; // use logsumexp?
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
    size_t n = tx.size()/q;// number of inputs

    assert(tx.size() >= 2 * q);
    assert(ty.size() >= 1 * q);
    assert(px.size() >= 2 * q);
    assert(py.size() >= 1 * q);

    std::fill(px.begin(), px.end(), 0.0);

    // Zeroth order (function value)
    std::vector<double> x0(n), x1(n), x2(n),xDivSum(n);
    for(size_t i = 0; i < n; ++i) {
        x0[i] = tx[i*q];  
    }

    auto max_iter = std::max_element(x0.begin(), x0.end());
    size_t max_idx = std::distance(x0.begin(), max_iter);
    double max_x = *max_iter;

    // Numerically stable log-sum-exp
    double sum_exp = 0.0;
    std::vector<double> expDiff(n);
    for (size_t i = 0; i < n; ++i) {
        expDiff[i] = std::exp(x0[i] - max_x);
        if (i == max_idx) continue; // skip max for stability
        sum_exp += expDiff[i];
    }
    double logsumexp = max_x + std::log1p(sum_exp);
    for (size_t i = 0; i < n; ++i) {
        xDivSum[i] = std::exp(x0[i]-logsumexp);
    }



// order 0
    if(order_up >= 0) {

    for(size_t i = 0; i < n; ++i)
        // ∂L/∂x0[i] += py[0] * ∂y0/∂x0[i] = py[0] * xDivSum[i]
        px[i*q] = py[0] * xDivSum[i];
    }

// order 1
    // S1 = sum x1[i]exp(x0[i]), S1div = S1/exp(maxx) S1divS = S1/S = S1div/sum(exp(xi-maxx))
    // S = sum exp(xi) = exp(maxx) sum(exp(xi - maxx))
    double S1div = 0.0, S1divS;
    if(order_up >= 1) {
        for(size_t i = 0; i < n; ++i) {
            S1div += expDiff[i] * tx[i*q+1];
        }
        S1divS = S1div/(sum_exp+1);
        for(size_t i = 0; i < n; ++i) {
            px[q*i] += py[1] * xDivSum[i] * (tx[i*q+1] - S1divS);
            px[q * i + 1] = py[1] * xDivSum[i];
        }
    }

    // not tested yet
    if(order_up >= 2) {
        for(size_t i=0; i<n; ++i) {
            // x2 direction
            px[i*q + 2] += py[2] * xDivSum[i];

            // accumulate x1 part and x0 part
            double contrib0 = 0.0;
            double contrib1 = 0.0;
            for(size_t j=0; j<n; ++j) {
                double H_ij = xDivSum[i] * ((i==j ? 1.0 : 0.0) - xDivSum[j]);
                // to x1[j]:
                contrib1 += (H_ij + xDivSum[j]*((j==i?1.0:0.0)-xDivSum[i])) * tx[j*q + 1];
                // to x0[i] via x2 term:
                contrib0 += H_ij * tx[i*q + 2];
                // to x0[j] via x1*x1 term:
                contrib0 -= xDivSum[j]* (tx[j*q + 1]*tx[i*q + 1]);
            }
            // x1 part
            px[i*q + 1] += py[2] * contrib1;
            // x0 part
            px[i*q + 0] += py[2] * contrib0;
        }

    }

//--- order 3 ---------------------------------------------------------------
    if(order_up >= 3) {
        // forward did:
        //   y3 = sum_i w[i] * x3[i]
        //      + 3 * sum_i w[i]*(x1[i]-μ1) * x2[i]
        //      + sum_i w[i]*(x1[i]-μ1)^3
        // where μ1 = sum_j w[j]*x1[j].
        //
        // reverse of that:
        //  ∂L/∂x3[i] += py[3] * w[i]
        //  ∂L/∂x2[i] += py[3] * 3 * w[i] * (x1[i] - μ1)
        //  ∂L/∂x1[i] += py[3] * [3 * w[i] * x2[i] + 3 * w[i] * (x1[i] - μ1)^2]
        //  ∂L/∂x0 via all w[i] dependencies — a bit messy, but you can
        //    reuse the derivative of w wrt x0: ∂w[i]/∂x0[k] = H_{ik}.

        // precompute μ1
        double mu1 = 0.0;
        for(size_t i=0; i<n; ++i) {
            mu1 += xDivSum[i] * tx[i*q + 1];
        }
        for(size_t i=0; i<n; ++i) {
            double dx1 = tx[i*q + 1] - mu1;
            // x3 part
            px[i*q + 3] += py[3] * xDivSum[i];
            // x2 part
            px[i*q + 2] += py[3] * 3.0 * xDivSum[i] * dx1;
            // x1 part
            px[i*q + 1] += py[3] * (3.0 * xDivSum[i] * tx[i*q + 2]
                                   + 3.0 * xDivSum[i] * dx1 * dx1);

            // now the contributions back to x0 via w[i] dependence
            // ∂w[i]/∂x0[k] = w[i]*((i==k?1.0:0.0) - w[k])
            for(size_t k=0; k<n; ++k) {
                double dwi = xDivSum[i] * ((i==k?1.0:0.0) - xDivSum[k]);
                double term = 0.0;
                // from x3 term:
                term += py[3] * tx[i*q + 3];
                // from x2 term:
                term += py[3] * 3.0 * dx1 * tx[i*q + 2];
                // from x1^3 term:
                term += py[3] * dx1 * dx1 * dx1;
                px[k*q + 0] += dwi * term;
            }
        }
    }

    return true;
}

atomic_logspace_add logspace_add_atomic("logspace_add_n");


template <class Type>
Type logspace_add_n(const CppAD::vector<Type>& x) {
    CppAD::vector<Type> yout(1);
    logspace_add_atomic(x, yout);
    return yout[0];
}

// Explicit instantiation for AD<double> and AD<AD<double>> if you want:
template CppAD::AD<double> logspace_add_n(const CppAD::vector<CppAD::AD<double>>&);
#endif
