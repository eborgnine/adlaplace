#ifndef ADLAPLACE_OMPAD_HPP
#define ADLAPLACE_OMPAD_HPP

#include <cstddef>

// Setup/teardown must run on OpenMP thread 0, outside any #pragma omp parallel
// region. Teardown flushes per-thread CppAD thread_alloc pools before
// parallel_setup(1) so a later CppadParallelScope (e.g. trace_hinv_t) is safe.

void cppad_parallel_setup(std::size_t num_threads);
void cppad_parallel_teardown();

struct CppadParallelScope {
  explicit CppadParallelScope(std::size_t num_threads) {
    cppad_parallel_setup(num_threads);
  }
  ~CppadParallelScope() { cppad_parallel_teardown(); }

  CppadParallelScope(const CppadParallelScope&) = delete;
  CppadParallelScope& operator=(const CppadParallelScope&) = delete;
};

#endif
