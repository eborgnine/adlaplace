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
  omp_set_dynamic(0);
  omp_set_num_threads(static_cast<int>(n));
}

static inline void cppad_parallel_setup(
  std::size_t num_threads
) {
  if (num_threads < 1) num_threads = 1;
  set_num_threads_wrapper(num_threads);
  CppAD::thread_alloc::parallel_setup(
    num_threads,
    &in_parallel_wrapper,
    &thread_num_wrapper
  );
  CppAD::parallel_ad<double>();
  CppAD::thread_alloc::hold_memory(false);
}

#endif // ADLAPLACE_OMPAD_HPP
