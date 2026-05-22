#include "adlaplace/ompad.hpp"

#include <omp.h>
#include <cppad/cppad.hpp>
#include <cppad/utility/thread_alloc.hpp>

namespace {

bool in_parallel_wrapper() {
  return omp_in_parallel() != 0;
}

std::size_t thread_num_wrapper() {
  return static_cast<std::size_t>(omp_get_thread_num());
}

void set_num_threads_wrapper(std::size_t n) {
  omp_set_dynamic(0);
  omp_set_num_threads(static_cast<int>(n));
}

}  // namespace

void cppad_parallel_setup(std::size_t num_threads) {
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
