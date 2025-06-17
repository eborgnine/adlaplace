// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <omp.h>

// Helper: for CppAD/Threading setup
bool in_parallel() { return omp_in_parallel() != 0; }
size_t thread_number() { return static_cast<size_t>(omp_get_thread_num()); }

// [[Rcpp::export]]
Rcpp::List parallel_jac_hess(Rcpp::NumericVector x, int num_threads = 2) {
    typedef CppAD::AD<double> ADdouble;
    typedef CppAD::vector<ADdouble> ADvec;
    typedef CppAD::vector<double> dvec;
    typedef CppAD::ADFun<double> ADFun;

    size_t nx = x.size();

    // 1. Build tape for function: y = prod(x)
    ADvec ax(nx), ay(1);
    for(size_t i = 0; i < nx; ++i)
        ax[i] = x[i];
    CppAD::Independent(ax);
    ay[0] = 0;
    for(size_t i = 0; i < nx; ++i)
        ay[0] += ax[i] * ax[i] * ax[i];
    ADFun fun(ax, ay);

    // 2. Threading setup
    omp_set_num_threads(num_threads);
    CppAD::thread_alloc::parallel_setup(num_threads, in_parallel, thread_number);
    CppAD::parallel_ad<double>();
    CppAD::thread_alloc::hold_memory(true);

    // 3. Value and Jacobian calculation
    dvec x0(nx);
    for(size_t i = 0; i < nx; ++i) x0[i] = x[i];
    std::vector<ADFun> fun_threads(num_threads);
    for(int i = 0; i < num_threads; ++i)
        fun_threads[i] = fun;

    double value = fun.Forward(0, x0)[0];

    Rcpp::NumericVector jacobian(nx);
    #pragma omp parallel for
    for(int j = 0; j < nx; ++j) {
        size_t tid = thread_number();
        dvec dx(nx, 0.0); dx[j] = 1.0;
        fun_threads[tid].Forward(0, x0); // clear Taylor memory
        dvec dy = fun_threads[tid].Forward(1, dx);
        jacobian[j] = dy[0];
    }

    // 4. Hessian calculation (parallelize by columns for safety)
    Rcpp::NumericMatrix hessian(nx, nx);
    #pragma omp parallel for
    for(int j = 0; j < nx; ++j) {
        size_t tid = thread_number();
        dvec w(1, 1.0);
        dvec u(nx, 0.0);
        u[j] = 1.0;
        fun_threads[tid].Forward(0, x0);
        fun_threads[tid].Forward(1, u);
        dvec ddw = fun_threads[tid].Reverse(2, w);
        for(size_t i = 0; i < nx; ++i)
            hessian(i, j) = ddw[2 * i + 1];
    }

    // 5. Cleanup
    CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
    CppAD::thread_alloc::hold_memory(false);
    for(int i = 1; i < num_threads; ++i)
        CppAD::thread_alloc::free_available(i);
    CppAD::thread_alloc::free_available(0);

    return Rcpp::List::create(
        Rcpp::Named("value") = value,
        Rcpp::Named("jacobian") = jacobian,
        Rcpp::Named("hessian") = hessian
    );
}


#ifdef UNDEF



library(Rcpp)
sourceCpp("hpol/src2/omp.cpp")
parallel_jac_hess(c(2, 3, 4), num_threads = 2)


#endif
