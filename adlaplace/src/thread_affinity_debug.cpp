#include "adlaplace/runtime/thread_affinity_debug.hpp"

#include <Rcpp.h>

#include <omp.h>

#include <sstream>
#include <string>
#include <vector>

#ifdef DEBUG

namespace {

struct ThreadMismatch {
  std::size_t shard = 0;
  std::size_t owner_thread = 0;
  std::size_t actual_thread = 0;
  std::string phase;
};

std::vector<ThreadMismatch>& mismatch_log() {
  static std::vector<ThreadMismatch> log;
  return log;
}

}  // namespace

bool adlaplace_debug_enabled() {
  return true;
}

bool adlaplace_shard_thread_ok(const GroupPack& pack) {
  if (!pack.owner_thread_assigned) {
    return true;
  }
  return static_cast<std::size_t>(omp_get_thread_num()) == pack.owner_thread;
}

void adlaplace_debug_note_grad_mismatch(double* grad_local, std::size_t n) {
  if (n == 0 || grad_local == nullptr) return;
#pragma omp critical(adlaplace_debug_grad)
  {
    grad_local[0] = adlaplace_debug::kThreadMismatchSentinel;
  }
}

void adlaplace_debug_note_trace_mismatch(double* trace_accum, std::size_t n) {
  if (n == 0 || trace_accum == nullptr) return;
#pragma omp critical(adlaplace_debug_trace)
  {
    trace_accum[0] = adlaplace_debug::kThreadMismatchSentinel;
  }
}

void adlaplace_debug_record_mismatch(
  std::size_t shard,
  std::size_t owner_thread,
  std::size_t actual_thread,
  const char* phase) {
  const char* phase_str = phase ? phase : "(unknown)";
#pragma omp critical(adlaplace_debug_mismatch_log)
  {
    mismatch_log().push_back(
      ThreadMismatch{shard, owner_thread, actual_thread, std::string(phase_str)}
    );
  }
}

void adlaplace_debug_raise_if_any(const char* context) {
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
  for (const ThreadMismatch& m : local) {
    os << "  shard " << m.shard
       << ": owner_thread=" << m.owner_thread
       << ", omp_get_thread_num()=" << m.actual_thread
       << ", phase=" << m.phase << "\n";
  }
  Rcpp::stop("%s", os.str().c_str());
}

void adlaplace_debug_print_load_banner() {
  static bool shown = false;
  if (shown) return;
  shown = true;
  Rcpp::Rcout << "adlaplace: DEBUG build (thread-affinity checks)\n";
}

#else

bool adlaplace_debug_enabled() {
  return false;
}

bool adlaplace_shard_thread_ok(const GroupPack&) {
  return true;
}

void adlaplace_debug_note_grad_mismatch(double*, std::size_t) {}

void adlaplace_debug_note_trace_mismatch(double*, std::size_t) {}

void adlaplace_debug_record_mismatch(
  std::size_t,
  std::size_t,
  std::size_t,
  const char*) {}

void adlaplace_debug_raise_if_any(const char*) {}

void adlaplace_debug_print_load_banner() {}

#endif
