#include "adlaplace/ompad.hpp"

#include <Rcpp.h>

#include <omp.h>
#include <cppad/cppad.hpp>
#include <cppad/utility/thread_alloc.hpp>

namespace {

std::size_t cppad_team_num_threads = 1;

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

void require_serial_main_thread(const char* phase) {
  if (omp_in_parallel() != 0) {
    Rcpp::stop("%s: must not run inside an OpenMP parallel region", phase);
  }
  if (omp_get_thread_num() != 0) {
    Rcpp::stop("%s: must run on OpenMP thread 0", phase);
  }
}

#ifdef DEBUG
void debug_teardown_flush(std::size_t n_flush) {
  Rcpp::Rcout << "cppad_parallel_teardown: flushing " << n_flush
              << " thread pools\n";
}
void debug_teardown_restored() {
  Rcpp::Rcout << "cppad_parallel_teardown: serial mode restored\n";
}
#else
void debug_teardown_flush(std::size_t) {}
void debug_teardown_restored() {}
#endif

}  // namespace

void cppad_parallel_setup(std::size_t num_threads) {
  require_serial_main_thread("cppad_parallel_setup");
  if (num_threads < 1) num_threads = 1;

  if (num_threads != cppad_team_num_threads) {
    cppad_parallel_teardown();
  }

  cppad_team_num_threads = num_threads;
  set_num_threads_wrapper(num_threads);

  if (num_threads == 1) {
    CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
    CppAD::thread_alloc::hold_memory(false);
  } else {
    CppAD::thread_alloc::parallel_setup(
      num_threads,
      &in_parallel_wrapper,
      &thread_num_wrapper
    );
    CppAD::thread_alloc::hold_memory(true);
  }
  CppAD::parallel_ad<double>();
}

void cppad_parallel_teardown() {
  require_serial_main_thread("cppad_parallel_teardown");

  const std::size_t n_flush = cppad_team_num_threads;
  if (n_flush > 1) {
    debug_teardown_flush(n_flush);
    set_num_threads_wrapper(n_flush);
#pragma omp parallel num_threads(static_cast<int>(n_flush))
    {
      CppAD::thread_alloc::free_available(
        static_cast<std::size_t>(omp_get_thread_num())
      );
    }
    for (std::size_t t = 0; t < n_flush; ++t) {
      CppAD::thread_alloc::free_available(t);
    }
    (void)CppAD::thread_alloc::free_all();
  }

  set_num_threads_wrapper(1);
  CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
  CppAD::thread_alloc::hold_memory(false);
  CppAD::parallel_ad<double>();
  cppad_team_num_threads = 1;
  debug_teardown_restored();
}
