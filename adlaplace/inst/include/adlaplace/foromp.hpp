// utilities for setting up omp threads

#ifndef ADLAPLACE_FOROMP_HPP
#define ADLAPLACE_FOROMP_HPP

#include <cstddef>

#include <cppad/utility/thread_alloc.hpp>
// If that header path doesn't exist in your setup, fall back to:
// #include <cppad/cppad.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

static inline bool in_parallel_wrapper() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

static inline std::size_t thread_num_wrapper() {
#ifdef _OPENMP
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

static inline std::size_t max_threads_wrapper() {
#ifdef _OPENMP
  return static_cast<std::size_t>(omp_get_max_threads());
#else
  return 1;
#endif
}

static inline void set_num_threads_wrapper(std::size_t n) {
#ifdef _OPENMP
  omp_set_num_threads(static_cast<int>(n));
#else
  (void)n;
#endif
}

// New: one-call setup for CppAD + OpenMP
static inline void cppad_parallel_setup(std::size_t nthreads) {
  set_num_threads_wrapper(nthreads);

  CppAD::thread_alloc::parallel_setup(
    nthreads,
    []() { return in_parallel_wrapper(); },
    []() { return thread_num_wrapper(); }
  );
}

#endif