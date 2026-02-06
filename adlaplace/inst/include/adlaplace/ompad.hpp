// functions to initiate omp threads compatible with cppad

#ifndef ADLAPLACE_FOROMP_HPP
#define ADLAPLACE_FOROMP_HPP

#include <cstddef>
#include <cppad/utility/thread_alloc.hpp>

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

static inline void set_num_threads_wrapper(std::size_t n) {
#ifdef _OPENMP
  omp_set_num_threads(static_cast<int>(n));
#else
  (void)n;
#endif
}

struct CppAD_OMP_Guard {
  std::size_t old_max_threads = 1;
  bool active = false;

  CppAD_OMP_Guard() = default;

  explicit CppAD_OMP_Guard(std::size_t nthreads) : active(true) {
#ifdef _OPENMP
    old_max_threads = static_cast<std::size_t>(omp_get_max_threads());
#endif

    set_num_threads_wrapper(nthreads);

#ifdef _OPENMP
    omp_set_max_active_levels(1); // no omp_set_nested()
#endif

    CppAD::thread_alloc::parallel_setup(
      nthreads,
      &in_parallel_wrapper,
      &thread_num_wrapper
    );
  }

  CppAD_OMP_Guard(const CppAD_OMP_Guard&) = delete;
  CppAD_OMP_Guard& operator=(const CppAD_OMP_Guard&) = delete;
  CppAD_OMP_Guard(CppAD_OMP_Guard&&) = delete;
  CppAD_OMP_Guard& operator=(CppAD_OMP_Guard&&) = delete;

  ~CppAD_OMP_Guard() {
    if (!active) return;

    CppAD::thread_alloc::parallel_setup(1, &in_parallel_wrapper, &thread_num_wrapper);
    CppAD::thread_alloc::free_all();

#ifdef _OPENMP
    omp_set_num_threads(static_cast<int>(old_max_threads));
    omp_set_max_active_levels(1);
#endif
  }
};

static inline CppAD_OMP_Guard cppad_parallel_setup(std::size_t nthreads) {
  return CppAD_OMP_Guard(nthreads);
}

#endif // ADLAPLACE_FOROMP_HPP