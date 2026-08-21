#include "adlaplace/omp_compat.hpp"
#include "adlaplace/ompad.hpp"
#include "adlaplace/runtime.hpp"

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <cppad/utility/thread_alloc.hpp>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::size_t cppad_team_num_threads = 1;

bool in_parallel_wrapper() { return omp_in_parallel() != 0; }

std::size_t thread_num_wrapper() {
  return static_cast<std::size_t>(omp_get_thread_num());
}

void set_num_threads_wrapper(std::size_t n) {
  omp_set_dynamic(0);
  omp_set_num_threads(static_cast<int>(n));
}

void require_serial_main_thread(const char *phase) {
  if (omp_in_parallel() != 0) {
    Rcpp::stop("%s: must not run inside an OpenMP parallel region", phase);
  }
  if (omp_get_thread_num() != 0) {
    Rcpp::stop("%s: must run on OpenMP thread 0", phase);
  }
}

// ADLAPLACE_HOLD_MEMORY=0 disables hold_memory(true) so return_memory deletes
// via the system allocator (thread-safe). Default is hold (1). Mitigation only;
// prefer draining shard buffers at parallel boundaries.
bool hold_memory_enabled() {
  const char *env = std::getenv("ADLAPLACE_HOLD_MEMORY");
  if (env == nullptr || env[0] == '\0') {
    return true;
  }
  return !(env[0] == '0' && env[1] == '\0');
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

} // namespace

void cppad_parallel_setup(std::size_t num_threads) {
  require_serial_main_thread("cppad_parallel_setup");
  if (num_threads < 1)
    num_threads = 1;
#ifndef _OPENMP
  num_threads = 1;
#endif

  if (num_threads != cppad_team_num_threads) {
    cppad_parallel_teardown();
  }

  cppad_team_num_threads = num_threads;
  set_num_threads_wrapper(num_threads);

  if (num_threads == 1) {
    CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
    CppAD::thread_alloc::hold_memory(false);
  } else {
    CppAD::thread_alloc::parallel_setup(num_threads, &in_parallel_wrapper,
                                        &thread_num_wrapper);
    CppAD::thread_alloc::hold_memory(hold_memory_enabled());
  }
  CppAD::parallel_ad<double>();
}

//' Whether this build was compiled with OpenMP support.
//'
//' @return \code{TRUE} if OpenMP was enabled at compile time, otherwise
//'   \code{FALSE}.
//' @export
// [[Rcpp::export(rng = false)]]
bool has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

//' Touch OpenMP from the main thread after \code{dyn.load}.
//'
//' Intended for fresh R processes (e.g. \code{R CMD check} vignette re-builds)
//' so Homebrew libomp TLS is initialized before the first multi-thread
//' CppAD session.
//'
//' @keywords internal
// [[Rcpp::export(rng = false)]]
void warm_openmp_runtime() {
#ifdef _OPENMP
  omp_set_dynamic(0);
  (void)omp_get_max_threads();
#pragma omp parallel num_threads(1)
  {
    (void)omp_get_thread_num();
  }
#endif
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
          static_cast<std::size_t>(omp_get_thread_num()));
    }
    for (std::size_t t = 0; t < n_flush; ++t) {
      CppAD::thread_alloc::free_available(t);
    }
    const bool all_freed = CppAD::thread_alloc::free_all();
#ifdef DEBUG
    if (!all_freed) {
      Rcpp::Rcout
          << "cppad_parallel_teardown: free_all() returned false "
             "(in-use blocks remain; check cross-thread thread_alloc frees)\n";
    }
#else
    (void)all_freed;
#endif
  }

  set_num_threads_wrapper(1);
  CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
  CppAD::thread_alloc::hold_memory(false);
  CppAD::parallel_ad<double>();
  cppad_team_num_threads = 1;
  debug_teardown_restored();
}

#ifdef DEBUG

namespace {

struct ThreadMismatch {
  std::size_t shard = 0;
  std::size_t owner_thread = 0;
  std::size_t actual_thread = 0;
  std::string phase;
};

std::vector<ThreadMismatch> &mismatch_log() {
  static std::vector<ThreadMismatch> log;
  return log;
}

} // namespace

bool adlaplace_debug_enabled() { return true; }

bool adlaplace_shard_thread_ok(const AdTape &pack) {
  if (!pack.owner_thread_assigned) {
    return true;
  }
  return static_cast<std::size_t>(omp_get_thread_num()) == pack.owner_thread;
}

void adlaplace_debug_note_grad_mismatch(double *grad_local, std::size_t n) {
  if (n == 0 || grad_local == nullptr)
    return;
#pragma omp critical(adlaplace_debug_grad)
  {
    grad_local[0] = adlaplace_debug::kThreadMismatchSentinel;
  }
}

void adlaplace_debug_note_trace_mismatch(double *trace_accum, std::size_t n) {
  if (n == 0 || trace_accum == nullptr)
    return;
#pragma omp critical(adlaplace_debug_trace)
  {
    trace_accum[0] = adlaplace_debug::kThreadMismatchSentinel;
  }
}

void adlaplace_debug_record_mismatch(std::size_t shard,
                                     std::size_t owner_thread,
                                     std::size_t actual_thread,
                                     const char *phase) {
  const char *phase_str = phase ? phase : "(unknown)";
#pragma omp critical(adlaplace_debug_mismatch_log)
  {
    mismatch_log().push_back(ThreadMismatch{shard, owner_thread, actual_thread,
                                            std::string(phase_str)});
  }
}

void adlaplace_debug_raise_if_any(const char *context) {
  std::vector<ThreadMismatch> local;
#pragma omp critical(adlaplace_debug_mismatch_log)
  {
    local.swap(mismatch_log());
  }
  if (local.empty()) {
    return;
  }

  std::ostringstream os;
  os << "OpenMP thread affinity mismatch";
  if (context && context[0] != '\0') {
    os << " (" << context << ")";
  }
  os << ":\n";
  for (const ThreadMismatch &m : local) {
    os << "  shard " << m.shard << ": owner_thread=" << m.owner_thread
       << ", omp_get_thread_num()=" << m.actual_thread << ", phase=" << m.phase
       << "\n";
  }
  Rcpp::stop("%s", os.str().c_str());
}

void adlaplace_debug_print_load_banner() {
  static bool shown = false;
  if (shown)
    return;
  shown = true;
  Rcpp::Rcout << "adlaplace: DEBUG build (thread-affinity checks)\n";
}

#else

bool adlaplace_debug_enabled() { return false; }

bool adlaplace_shard_thread_ok(const AdTape &) { return true; }

void adlaplace_debug_note_grad_mismatch(double *, std::size_t) {}

void adlaplace_debug_note_trace_mismatch(double *, std::size_t) {}

void adlaplace_debug_record_mismatch(std::size_t, std::size_t, std::size_t,
                                     const char *) {}

void adlaplace_debug_raise_if_any(const char *) {}

void adlaplace_debug_print_load_banner() {}

#endif
