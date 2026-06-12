#ifndef ADLAPLACE_OMPAD_HPP
#define ADLAPLACE_OMPAD_HPP

#include <cstddef>

// Setup/teardown must run on OpenMP thread 0, outside any #pragma omp parallel
// region. Tape build stays in default sequential CppAD mode; CppadParallelScope(N)
// at eval time (inner_opt, trace_hinv_t, etc.) owns parallel setup/teardown.
// cppad_parallel_setup() tears down first when the thread count changes.

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
