// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <omp.h>

// Helper: for CppAD/Threading setup
bool in_parallel() { return omp_in_parallel() != 0; }
size_t thread_number() { return static_cast<size_t>(omp_get_thread_num()); }

// [[Rcpp::export]]
Rcpp::NumericVector parallel_jacobian(Rcpp::NumericVector x, int num_threads = 2) {
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
    ay[0] = ax[0];
    for(size_t i = 1; i < nx; ++i)
        ay[0] *= ax[i];
    ADFun fun(ax, ay);

    // 2. Threading setup
    omp_set_num_threads(num_threads);
    CppAD::thread_alloc::parallel_setup(num_threads, in_parallel, thread_number);
    CppAD::parallel_ad<double>();
    CppAD::thread_alloc::hold_memory(true);

    // 3. Parallel Jacobian calculation
    Rcpp::NumericVector jacobian(nx);
    dvec x0(nx);
    for(size_t i = 0; i < nx; ++i) x0[i] = x[i];

    std::vector<ADFun> fun_threads(num_threads);
    for(int i = 0; i < num_threads; ++i)
        fun_threads[i] = fun;

    #pragma omp parallel for
    for(int j = 0; j < nx; ++j) {
        size_t tid = thread_number();
        dvec dx(nx, 0.0); dx[j] = 1.0;
        fun_threads[tid].Forward(0, x0); // clear Taylor memory
        dvec dy = fun_threads[tid].Forward(1, dx);
        jacobian[j] = dy[0];
    }

    // 4. Cleanup
    CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
    CppAD::thread_alloc::hold_memory(false);
    for(int i = 1; i < num_threads; ++i)
        CppAD::thread_alloc::free_available(i);
    CppAD::thread_alloc::free_available(0);

    return jacobian;
}

#ifdef UNDEF

Sys.setenv("PKG_CPPFLAGS" = "-I/opt/homebrew/include")

library(Rcpp)
sourceCpp("hpol/src2/omp.cpp")

x <- c(2, 3, 4)
parallel_jacobian(x, num_threads = 4)
# Should return c(3*4, 2*4, 2*3) = c(12, 8, 6)

#endif
