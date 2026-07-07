#ifndef ADLAPLACE_TRACE_HINV_T_RUNTIME_HPP
#define ADLAPLACE_TRACE_HINV_T_RUNTIME_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <omp.h>
#include <vector>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/ompad.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"
#include "adlaplace/runtime/thread_affinity_debug.hpp"
#include "adlaplace/runtime/thread_groups.hpp"

namespace adlaplace_trace {

inline void reset_shard_adfun_taylor(GroupPack& gp) {
  gp.fun.capacity_order(0);
}

inline void trace_hinv_t_reset_shards(ad_fun& backend) {
  for (std::size_t s = 0; s < backend.fun.size(); ++s) {
    GroupPack& gp = *pack_ctx(shard_handle(&backend, s)->ctx);
    reset_shard_adfun_taylor(gp);
    gp.direction.clear();
    gp.direction_zeros.clear();
    gp.wthree.clear();
  }
}

// Assumes an active CppadParallelScope; does not open a nested scope.
inline std::vector<double> trace_hinv_t_parallel(
  ad_fun& backend,
  const std::vector<double>& x,
  const std::vector<int>& LinvPt_p,
  const std::vector<int>& LinvPt_i,
  const std::vector<double>& LinvPt_x,
  const std::size_t LinvPt_ncol,
  const std::vector<int>& LinvPtColumns_p,
  const std::vector<int>& LinvPtColumns_i,
  const bool verbose) {
    
  const std::size_t Nparams = pack_ctx(backend.fun[0]->ctx)->x.size();
  if (x.size() != Nparams) {
    Rcpp::stop(
      "x has length %d but expected Nparams=%d",
      static_cast<int>(x.size()),
      static_cast<int>(Nparams)
    );
  }

  const std::size_t LinvPt_p_len = LinvPt_p.size();
  const std::size_t LinvPt_i_len = LinvPt_i.size();
  const std::size_t LinvPt_x_len = LinvPt_x.size();
  const std::size_t LinvPtColumns_p_len = LinvPtColumns_p.size();
  const std::size_t LinvPtColumns_i_len = LinvPtColumns_i.size();

  std::vector<double> trace_accum(Nparams, 0.0);
  if (backend.fun.empty()) {
    return trace_accum;
  }

  adlaplace_require_owner_threads_assigned(backend);
  const std::vector<std::vector<std::size_t>> thread_groups =
    thread_groups_from_backend(backend);
  const int num_threads = static_cast<int>(thread_groups.size());

  if (verbose) {
    Rcpp::Rcout << "trace_hinv_t: threads = " << num_threads
                << ", shards = " << backend.fun.size()
                << ", params = " << Nparams << "\n";
  }

#pragma omp parallel num_threads(num_threads)
  {
    const int tid = omp_get_thread_num();
    std::vector<double> trace_local(Nparams, 0.0);
    const std::vector<std::size_t>& shard_group =
      thread_groups[static_cast<std::size_t>(tid)];
    for (std::size_t s : shard_group) {
      adlaplace_adpack_handle* h = shard_handle(&backend, s);
      h->api->assign_memory(h->ctx);
      h->api->trace_hinv_t(
        h->ctx,
        x.data(),
        LinvPt_p.data(),
        LinvPt_i.data(),
        LinvPt_x.data(),
        LinvPt_ncol,
        LinvPt_p_len,
        LinvPt_i_len,
        LinvPt_x_len,
        LinvPtColumns_p.data(),
        LinvPtColumns_i.data(),
        LinvPtColumns_p_len,
        LinvPtColumns_i_len,
        trace_local.data()
      );
    }
#pragma omp critical(trace_hinv_merge)
    {
      for (std::size_t d = 0; d < Nparams; ++d) {
        trace_accum[d] += trace_local[d];
      }
    }
  }

  if (verbose) {
    Rcpp::Rcout << "trace_hinv_t: finished\n";
  }

  return trace_accum;
}

inline std::vector<double> trace_hinv_t_impl(
  ad_fun& backend,
  const std::vector<double>& x,
  const std::vector<int>& LinvPt_p,
  const std::vector<int>& LinvPt_i,
  const std::vector<double>& LinvPt_x,
  const std::size_t LinvPt_ncol,
  const std::vector<int>& LinvPtColumns_p,
  const std::vector<int>& LinvPtColumns_i,
  const bool verbose) {
  if (backend.fun.empty()) {
    return std::vector<double>();
  }

  adlaplace_require_owner_threads_assigned(backend);
  const std::vector<std::vector<std::size_t>> thread_groups =
    thread_groups_from_backend(backend);
  const int num_threads = static_cast<int>(thread_groups.size());

  if (verbose && adlaplace_debug_enabled()) {
    Rcpp::Rcout << "trace_hinv_t: entering CppAD parallel session ("
                << num_threads << " threads)\n";
  }

  trace_hinv_t_reset_shards(backend);

  std::vector<double> trace_accum;
  {
    CppadParallelScope parallel_scope(static_cast<std::size_t>(num_threads));
    trace_accum = trace_hinv_t_parallel(
      backend,
      x,
      LinvPt_p,
      LinvPt_i,
      LinvPt_x,
      LinvPt_ncol,
      LinvPtColumns_p,
      LinvPtColumns_i,
      verbose
    );
    if (verbose && adlaplace_debug_enabled()) {
      Rcpp::Rcout << "trace_hinv_t: shard eval done\n";
    }
  }

  if (verbose && adlaplace_debug_enabled()) {
    Rcpp::Rcout << "trace_hinv_t: parallel block ended\n";
  }
  adlaplace_debug_raise_if_any("trace_hinv_t after parallel block");

  return trace_accum;
}

}  // namespace adlaplace_trace

#endif
