#ifndef ADLAPLACE_FOROMP_HPP
#define ADLAPLACE_FOROMP_HPP

#include <cstddef>
#include <cppad/utility/thread_alloc.hpp>

#ifdef _OPENMP
  #include <omp.h>
#endif

// ---- wrappers as plain function pointers (CppAD wants function pointers) ----
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

// ---- guard object: sets up in ctor, resets in dtor ----
struct CppAD_OMP_Guard {
  std::size_t old_max_threads = 1;
  bool active = false;

  CppAD_OMP_Guard() = default;

  explicit CppAD_OMP_Guard(std::size_t nthreads)
  : active(true)
  {
#ifdef _OPENMP
    // Save current setting so we can restore it
    old_max_threads = static_cast<std::size_t>(omp_get_max_threads());
#endif

    // Set OpenMP to requested threads
    set_num_threads_wrapper(nthreads);

#ifdef _OPENMP
    // Avoid nested parallelism surprises
    omp_set_nested(0);
    omp_set_max_active_levels(1);
#endif

    // Tell CppAD about the OpenMP threading model
    CppAD::thread_alloc::parallel_setup(
      nthreads,
      &in_parallel_wrapper,
      &thread_num_wrapper
    );
  }

  // no copy
  CppAD_OMP_Guard(const CppAD_OMP_Guard&) = delete;
  CppAD_OMP_Guard& operator=(const CppAD_OMP_Guard&) = delete;

  // movable
  CppAD_OMP_Guard(CppAD_OMP_Guard&& other) noexcept
  : old_max_threads(other.old_max_threads), active(other.active)
  {
    other.active = false;
  }
  CppAD_OMP_Guard& operator=(CppAD_OMP_Guard&& other) noexcept {
    if (this != &other) {
      // finish current guard first
      if (active) this->~CppAD_OMP_Guard();
      old_max_threads = other.old_max_threads;
      active = other.active;
      other.active = false;
    }
    return *this;
  }

  ~CppAD_OMP_Guard() {
    if (!active) return;

    // Put CppAD back into single-thread mode and free thread-local alloc
    CppAD::thread_alloc::parallel_setup(1, &in_parallel_wrapper, &thread_num_wrapper);
    CppAD::thread_alloc::free_all();

#ifdef _OPENMP
    // Restore OpenMP thread setting
    omp_set_num_threads(static_cast<int>(old_max_threads));
    // Keep nesting off (generally safest inside R packages)
    omp_set_nested(0);
    omp_set_max_active_levels(1);
#endif
  }
};

// Factory function: returns guard (RAII)
static inline CppAD_OMP_Guard cppad_parallel_setup(std::size_t nthreads) {
  return CppAD_OMP_Guard(nthreads);
}

#endif // ADLAPLACE_FOROMP_HPP

