#ifndef ADLAPLACE_RUNTIME_HPP
#define ADLAPLACE_RUNTIME_HPP

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "adlaplace/backend.hpp"
#include "adlaplace/rviews.hpp"

ad_fun* clone_ad_fun(const ad_fun* src);

hessian_template hessian_template_from_dgc(
  const DgCView& tpl,
  const char* label);

void ad_fun_attach_hessians_from_list(
  ad_fun& shards,
  const Rcpp::List& ad_fun);

adlaplace_shard* shard_handle(ad_fun* backend, size_t shard);

SEXP ad_fun_handle_sexp(const Rcpp::List& ad_fun_list);

ad_fun* ad_fun_from_list(const Rcpp::List& ad_fun_list);

ad_fun* ad_fun_from_handle(SEXP handle);

ad_fun* resolve_ad_fun_eval(SEXP ad_fun_ptr);

ad_fun* resolve_ad_fun_laplace(const Rcpp::S4& ad_fun_s4);

std::vector<size_t> resolve_shard_indices(
  size_t n_shards,
  const Rcpp::IntegerVector& shards);

Rcpp::List sparsity_shard_from_handle(adlaplace_shard* shard);

inline void adlaplace_require_owner_threads_assigned(const ad_fun& backend) {
  for (std::size_t s = 0; s < backend.fun.size(); ++s) {
    adlaplace_shard* shard = shard_handle(const_cast<ad_fun*>(&backend), s);
    if (!shard->pack.owner_thread_assigned) {
      Rcpp::stop(
        "OpenMP thread assignment missing; call ad_fun(..., num_threads=) before inner_opt or trace_hinv_t"
      );
    }
  }
}

inline std::vector<std::vector<std::size_t>> thread_groups_from_backend(ad_fun& backend) {
  const std::size_t n = backend.fun.size();
  std::vector<std::size_t> owners(n);
  std::size_t max_t = 0;
  for (std::size_t s = 0; s < n; ++s) {
    const std::size_t t = backend.fun[s]->pack.owner_thread;
    owners[s] = t;
    if (t > max_t) max_t = t;
  }
  std::vector<std::vector<std::size_t>> groups(max_t + 1);
  for (std::size_t s = 0; s < n; ++s) {
    groups[owners[s]].push_back(s);
  }
  return groups;
}

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
