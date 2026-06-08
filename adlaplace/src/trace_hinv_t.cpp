#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <vector>

#include <omp.h>

#include "adlaplace/ompad.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"
#include "adlaplace/runtime/thread_affinity_debug.hpp"
#include "adlaplace/runtime/thread_groups.hpp"

namespace {

static void reset_shard_adfun_taylor(GroupPack& gp) {
  gp.fun.capacity_order(0);
}

int eval_shard_assign_memory(ad_fun& backend, std::size_t s) {
  adlaplace_adpack_handle* h = shard_handle(&backend, s);
  GroupPack& pack = *pack_ctx(h->ctx);
  if (!adlaplace_shard_thread_ok(pack)) {
    adlaplace_debug_record_mismatch(
      s,
      pack.owner_thread,
      static_cast<std::size_t>(omp_get_thread_num()),
      "trace_hinv_t assign_memory"
    );
    return 0;
  }
  if (!h->api->assign_memory) {
    return 18;
  }
  return h->api->assign_memory(h->ctx);
}

int eval_shard_trace_hinv_t(
  ad_fun& backend,
  std::size_t s,
  const double* x,
  const int* LinvPt_p,
  const int* LinvPt_i,
  const double* LinvPt_x,
  std::size_t LinvPt_ncol,
  std::size_t LinvPt_p_len,
  std::size_t LinvPt_i_len,
  std::size_t LinvPt_x_len,
  const int* LinvPtColumns_p,
  const int* LinvPtColumns_i,
  std::size_t LinvPtColumns_p_len,
  std::size_t LinvPtColumns_i_len,
  double* trace_accum) {
  adlaplace_adpack_handle* h = shard_handle(&backend, s);
  GroupPack& pack = *pack_ctx(h->ctx);
  if (!adlaplace_shard_thread_ok(pack)) {
    adlaplace_debug_record_mismatch(
      s,
      pack.owner_thread,
      static_cast<std::size_t>(omp_get_thread_num()),
      "trace_hinv_t"
    );
    adlaplace_debug_note_trace_mismatch(trace_accum, pack.x.size());
    return 0;
  }
  if (!h->api->trace_hinv_t) {
    return 1;
  }
  return h->api->trace_hinv_t(
    h->ctx,
    x,
    LinvPt_p,
    LinvPt_i,
    LinvPt_x,
    LinvPt_ncol,
    LinvPt_p_len,
    LinvPt_i_len,
    LinvPt_x_len,
    LinvPtColumns_p,
    LinvPtColumns_i,
    LinvPtColumns_p_len,
    LinvPtColumns_i_len,
    trace_accum
  );
}

std::vector<double> trace_hinv_t_impl(
  ad_fun& backend,
  const std::vector<double>& x,
  const std::vector<int>& LinvPt_p,
  const std::vector<int>& LinvPt_i,
  const std::vector<double>& LinvPt_x,
  const std::size_t LinvPt_ncol,
  const std::vector<int>& LinvPtColumns_p,
  const std::vector<int>& LinvPtColumns_i,
  bool verbose) {
  const std::size_t Nparams = pack_ctx(backend.fun[0]->ctx)->fun.Domain();
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
  if (verbose && adlaplace_debug_enabled()) {
    Rcpp::Rcout << "trace_hinv_t: DEBUG build (thread-affinity checks active)\n";
  }

  int rc_error = 0;
  int rc_group = -1;

  if (verbose && adlaplace_debug_enabled()) {
    Rcpp::Rcout << "trace_hinv_t: entering CppAD parallel session ("
                << num_threads << " threads)\n";
  }
  for (std::size_t s = 0; s < backend.fun.size(); ++s) {
    GroupPack& gp = *pack_ctx(shard_handle(&backend, s)->ctx);
    reset_shard_adfun_taylor(gp);
    gp.direction.clear();
    gp.direction_zeros.clear();
    gp.wthree.clear();
  }
  {
    CppadParallelScope parallel_scope(static_cast<std::size_t>(num_threads));
#pragma omp parallel num_threads(num_threads)
    {
      const int tid = omp_get_thread_num();
      if (tid < 0 || tid >= num_threads) {
#pragma omp critical
        {
          if (rc_error == 0) {
            rc_error = -1;
            rc_group = tid;
          }
        }
      } else {
        std::vector<double> trace_local(Nparams, 0.0);
        const std::vector<std::size_t>& shard_group =
          thread_groups[static_cast<std::size_t>(tid)];
        for (std::size_t s : shard_group) {
          const int rc_assign = eval_shard_assign_memory(backend, s);
          if (rc_assign != 0) {
#pragma omp critical
            {
              if (rc_error == 0) {
                rc_error = rc_assign;
                rc_group = static_cast<int>(s);
              }
            }
            continue;
          }
          const int rc = eval_shard_trace_hinv_t(
            backend,
            s,
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
          if (rc != 0) {
#pragma omp critical
            {
              if (rc_error == 0) {
                rc_error = rc;
                rc_group = static_cast<int>(s);
              }
            }
          }
        }
#pragma omp critical(trace_hinv_merge)
        {
          for (std::size_t d = 0; d < Nparams; ++d) {
            trace_accum[d] += trace_local[d];
          }
        }
      }
    }
    adlaplace_debug_raise_if_any("trace_hinv_t shard eval");
    if (verbose && adlaplace_debug_enabled()) {
      Rcpp::Rcout << "trace_hinv_t: shard eval done\n";
    }
  }

  if (verbose && adlaplace_debug_enabled()) {
    Rcpp::Rcout << "trace_hinv_t: parallel block ended\n";
  }
  adlaplace_debug_raise_if_any("trace_hinv_t after parallel block");

  if (rc_error != 0) {
    Rcpp::stop(
      "backend api->trace_hinv_t failed for shard %d (code %d: %s)",
      rc_group,
      rc_error,
      trace_hinv_t_strerror(rc_error)
    );
  }

  if (verbose) {
    Rcpp::Rcout << "trace_hinv_t: finished\n";
  }

  return trace_accum;
}

}  // namespace

//' @rdname adlaplace_cpp
//' @param verbose Logical; if \code{TRUE}, print threads, shards, and parameter sizes.
// [[Rcpp::export]]
Rcpp::NumericVector trace_hinv_t(
  SEXP ad_fun_ptr,
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  bool verbose = false
) {
  ad_fun* backend = resolve_ad_fun_eval(ad_fun_ptr);

  Rcpp::IntegerVector LinvPt_p = LinvPt.slot("p");
  Rcpp::IntegerVector LinvPt_i = LinvPt.slot("i");
  Rcpp::NumericVector LinvPt_x = LinvPt.slot("x");
  Rcpp::IntegerVector LinvPt_Dim = LinvPt.slot("Dim");

  Rcpp::IntegerVector LinvPtColumns_p = LinvPtColumns.slot("p");
  Rcpp::IntegerVector LinvPtColumns_i = LinvPtColumns.slot("i");

  const std::size_t LinvPt_ncol = static_cast<std::size_t>(LinvPt_Dim[1]);
  const std::vector<double> x_vec(x.begin(), x.end());
  const std::vector<int> LinvPt_p_vec(LinvPt_p.begin(), LinvPt_p.end());
  const std::vector<int> LinvPt_i_vec(LinvPt_i.begin(), LinvPt_i.end());
  const std::vector<double> LinvPt_x_vec(LinvPt_x.begin(), LinvPt_x.end());
  const std::vector<int> LinvPtColumns_p_vec(LinvPtColumns_p.begin(), LinvPtColumns_p.end());
  const std::vector<int> LinvPtColumns_i_vec(LinvPtColumns_i.begin(), LinvPtColumns_i.end());

  const std::vector<double> trace_accum = trace_hinv_t_impl(
    *backend,
    x_vec,
    LinvPt_p_vec,
    LinvPt_i_vec,
    LinvPt_x_vec,
    LinvPt_ncol,
    LinvPtColumns_p_vec,
    LinvPtColumns_i_vec,
    verbose
  );

  return Rcpp::NumericVector(trace_accum.begin(), trace_accum.end());
}
