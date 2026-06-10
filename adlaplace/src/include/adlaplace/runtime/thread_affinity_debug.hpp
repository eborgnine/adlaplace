#ifndef ADLAPLACE_RUNTIME_THREAD_AFFINITY_DEBUG_HPP
#define ADLAPLACE_RUNTIME_THREAD_AFFINITY_DEBUG_HPP

#include <cmath>
#include <cstddef>
#include <limits>

#include "adlaplace/api/backend.hpp"

namespace adlaplace_debug {

inline constexpr double kThreadMismatchSentinel =
  std::numeric_limits<double>::quiet_NaN();

}  // namespace adlaplace_debug

bool adlaplace_debug_enabled();

bool adlaplace_shard_thread_ok(const GroupPack& pack);

void adlaplace_debug_note_grad_mismatch(double* grad_local, std::size_t n);

void adlaplace_debug_note_trace_mismatch(double* trace_accum, std::size_t n);

void adlaplace_debug_record_mismatch(
  std::size_t shard,
  std::size_t owner_thread,
  std::size_t actual_thread,
  const char* phase);

void adlaplace_debug_raise_if_any(const char* context);

void adlaplace_debug_print_load_banner();

#endif
