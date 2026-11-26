#ifndef ADLAPLACE_PARALLEL_BUILD_HPP
#define ADLAPLACE_PARALLEL_BUILD_HPP

// Header-only CppAD parallel-mode guard for tape *build* (taping +
// sparsity calibration) in whichever DSO includes it.
//
// Why header-only / per-DSO: build_ad_fun_obs lives in ad_pack.hpp and is
// compiled into adlaplace.so and into every backend .so. Backends do not
// link adlaplace.so, and CppAD thread_alloc / AD statics are per shared
// library, so the CppadParallelScope in adlaplace/src/threads.cpp (which
// owns *evaluation* in adlaplace.so) cannot set up a backend's CppAD.
// This guard performs the equivalent setup/teardown using the including
// DSO's own statics, mirroring cppad_parallel_setup/cppad_parallel_teardown.
//
// Must be constructed and destroyed on the main thread, outside any
// OpenMP parallel region and while no tape is recording.

#include <cppad/cppad.hpp>
#include <cppad/utility/thread_alloc.hpp>

#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace adlaplace_parallel_build {

// No serial omp_* stubs here. adlaplace/src/include/adlaplace/omp_compat.hpp
// already provides those for this package, and a second definition fails
// the no-OpenMP builds (macOS, wasm) that include both headers. These
// callbacks are only passed to CppAD when OpenMP is on.
inline bool in_parallel_callback() {
#ifdef _OPENMP
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

inline std::size_t thread_num_callback() {
#ifdef _OPENMP
  return static_cast<std::size_t>(omp_get_thread_num());
#else
  return 0;
#endif
}

// Requested thread count clamped to a usable team size: at least 1, at
// most the number of groups, and 1 when OpenMP is unavailable.
inline std::size_t effective_threads(int requested, std::size_t n_groups) {
#ifndef _OPENMP
  (void)requested;
  (void)n_groups;
  return 1;
#else
  if (requested < 1 || n_groups < 2) {
    return 1;
  }
  std::size_t n = static_cast<std::size_t>(requested);
  if (n > n_groups) {
    n = n_groups;
  }
  // Cap by hardware concurrency, NOT omp_get_max_threads(): teardown
  // resets the OpenMP thread limit to 1, which would otherwise
  // permanently serialise every later build in the process.
  const int num_procs = omp_get_num_procs();
  if (num_procs > 0 && n > static_cast<std::size_t>(num_procs)) {
    n = static_cast<std::size_t>(num_procs);
  }
  return n;
#endif
}

inline void setup(std::size_t num_threads) {
  if (num_threads <= 1) {
    return;
  }
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(static_cast<int>(num_threads));
  CppAD::thread_alloc::parallel_setup(num_threads, &in_parallel_callback,
                                      &thread_num_callback);
  CppAD::thread_alloc::hold_memory(true);
  CppAD::parallel_ad<double>();
#endif
}

inline void teardown(std::size_t num_threads) {
  if (num_threads <= 1) {
    return;
  }
#ifdef _OPENMP
  // Flush per-thread free lists for every team member while the tapes
  // built on those threads are still alive (in-use blocks are left
  // alone), then restore the default sequential CppAD mode, matching
  // cppad_parallel_teardown() in adlaplace/src/threads.cpp.
  omp_set_num_threads(static_cast<int>(num_threads));
#pragma omp parallel num_threads(static_cast<int>(num_threads))
  {
    CppAD::thread_alloc::free_available(
        static_cast<std::size_t>(omp_get_thread_num()));
  }
  for (std::size_t t = 0; t < num_threads; ++t) {
    CppAD::thread_alloc::free_available(t);
  }
  (void)CppAD::thread_alloc::free_all();
  omp_set_num_threads(1);
  CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
  CppAD::thread_alloc::hold_memory(false);
  CppAD::parallel_ad<double>();
#endif
}

// RAII wrapper: parallel CppAD mode for the lifetime of the guard.
// Inactive (a no-op) for a single thread.
struct Guard {
  explicit Guard(std::size_t num_threads)
      : num_threads_(num_threads), active_(num_threads > 1) {
    if (active_) {
      setup(num_threads_);
    }
  }
  ~Guard() {
    if (active_) {
      teardown(num_threads_);
    }
  }
  Guard(const Guard &) = delete;
  Guard &operator=(const Guard &) = delete;

  bool active() const { return active_; }
  std::size_t num_threads() const { return num_threads_; }

private:
  std::size_t num_threads_;
  bool active_;
};

} // namespace adlaplace_parallel_build

#endif
