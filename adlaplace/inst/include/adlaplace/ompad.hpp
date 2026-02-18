#ifndef ADLAPLACE_OMPAD_HPP
#define ADLAPLACE_OMPAD_HPP

#include <cstddef>
#include <omp.h>
#include <cppad/utility/thread_alloc.hpp>

static inline bool in_parallel_wrapper() { return omp_in_parallel() != 0; }
static inline std::size_t thread_num_wrapper() {
  return static_cast<std::size_t>(omp_get_thread_num());
}

static inline void set_num_threads_wrapper(std::size_t n) {
  omp_set_num_threads(static_cast<int>(n));
}

// Compatibility shim: keep the existing API with minimal one-way setup
// (no teardown/restoration behavior).
struct CppAD_OMP_Guard {
  CppAD_OMP_Guard() = default;
  CppAD_OMP_Guard(std::size_t cppad_threads, std::size_t omp_threads) {
    if (cppad_threads < 1) cppad_threads = 1;
    if (omp_threads < 1) omp_threads = 1;
    if (cppad_threads < omp_threads) cppad_threads = omp_threads;
    set_num_threads_wrapper(omp_threads);
    CppAD::thread_alloc::parallel_setup(
      cppad_threads,
      &in_parallel_wrapper,
      &thread_num_wrapper
    );
    CppAD::parallel_ad<double>();
    CppAD::thread_alloc::hold_memory(true);
  }
};

static inline CppAD_OMP_Guard cppad_parallel_setup(
  std::size_t cppad_threads,
  std::size_t omp_threads
) {
  return CppAD_OMP_Guard(cppad_threads, omp_threads);
}

#endif // ADLAPLACE_OMPAD_HPP
