#ifndef ADLAPLACE_FOROMP_HPP
#define ADLAPLACE_FOROMP_HPP

#include <cstddef>  // for size_t

#ifdef _OPENMP
  #include <omp.h>

  static inline bool in_parallel_wrapper() {
    return omp_in_parallel() != 0;
  }

  static inline std::size_t thread_num_wrapper() {
    return static_cast<std::size_t>(omp_get_thread_num());
  }

  static inline std::size_t max_threads_wrapper() {
    return static_cast<std::size_t>(omp_get_max_threads());
  }

  static inline void set_num_threads_wrapper(std::size_t n) {
    // omp_set_num_threads takes int
    omp_set_num_threads(static_cast<int>(n));
  }


#else
  // OpenMP not enabled: provide single-thread fallbacks
  static inline bool in_parallel_wrapper() { return false; }
  static inline std::size_t thread_num_wrapper() { return 0; }
  static inline std::size_t max_threads_wrapper() { return 1; }
  static inline void set_num_threads_wrapper(std::size_t) { /* no-op */ }
#endif

#endif // ADLAPLACE_FOROMP_HPP