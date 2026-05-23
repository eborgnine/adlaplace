#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <vector>

#include "adlaplace/ompad.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"

namespace {

Rcpp::IntegerVector shards_vector(
  const Rcpp::Nullable<Rcpp::IntegerVector>& shards) {
  if (shards.isNull()) {
    return Rcpp::IntegerVector();
  }
  return shards.get();
}

}  // namespace

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
  SEXP ad_fun_ptr,
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  const int num_threads,
  Rcpp::Nullable<Rcpp::IntegerVector> shards = R_NilValue
) {
  ad_fun* backend = resolve_ad_fun_eval(ad_fun_ptr);
  const size_t n_shards = backend->fun.size();
  const size_t Nparams = pack_ctx(backend->fun[0]->ctx)->x.size();
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  Rcpp::IntegerVector LinvPt_p = LinvPt.slot("p");
  Rcpp::IntegerVector LinvPt_i = LinvPt.slot("i");
  Rcpp::NumericVector LinvPt_x = LinvPt.slot("x");
  Rcpp::IntegerVector LinvPt_Dim = LinvPt.slot("Dim");

  Rcpp::IntegerVector LinvPtColumns_p = LinvPtColumns.slot("p");
  Rcpp::IntegerVector LinvPtColumns_i = LinvPtColumns.slot("i");

  const size_t LinvPt_ncol = static_cast<size_t>(LinvPt_Dim[1]);
  const size_t LinvPt_p_len = static_cast<size_t>(LinvPt_p.size());
  const size_t LinvPt_i_len = static_cast<size_t>(LinvPt_i.size());
  const size_t LinvPt_x_len = static_cast<size_t>(LinvPt_x.size());
  const size_t LinvPtColumns_p_len = static_cast<size_t>(LinvPtColumns_p.size());
  const size_t LinvPtColumns_i_len = static_cast<size_t>(LinvPtColumns_i.size());

  std::vector<double> trace_accum(Nparams, 0.0);
  const std::vector<size_t> shard_idx =
    resolve_shard_indices(n_shards, shards_vector(shards));

  int rc_error = 0;
  int rc_group = -1;
  const int n_shards_loop = static_cast<int>(shard_idx.size());

  if (n_shards_loop > 0) {
    cppad_parallel_setup(static_cast<std::size_t>(num_threads));

#pragma omp parallel num_threads(num_threads)
    {
      std::vector<double> trace_local(Nparams, 0.0);
      int rc_local = 0;
      int rc_group_local = -1;

#pragma omp for schedule(static,1)
      for (int sidx = 0; sidx < n_shards_loop; ++sidx) {
        const size_t s = shard_idx[static_cast<std::size_t>(sidx)];
        adlaplace_adpack_handle* h = shard_handle(backend, s);
        if (!h->api->trace_hinv_t) {
          rc_local = 1;
          rc_group_local = static_cast<int>(s);
          continue;
        }
        const int rc = h->api->trace_hinv_t(
          h->ctx,
          x.begin(),
          LinvPt_p.begin(),
          LinvPt_i.begin(),
          LinvPt_x.begin(),
          LinvPt_ncol,
          LinvPt_p_len,
          LinvPt_i_len,
          LinvPt_x_len,
          LinvPtColumns_p.begin(),
          LinvPtColumns_i.begin(),
          LinvPtColumns_p_len,
          LinvPtColumns_i_len,
          trace_local.data()
        );
        if (rc != 0 && rc_local == 0) {
          rc_local = rc;
          rc_group_local = static_cast<int>(s);
        }
      }

#pragma omp critical
      {
        if (rc_error == 0 && rc_local != 0) {
          rc_error = rc_local;
          rc_group = rc_group_local;
        }
        for (size_t d = 0; d < Nparams; ++d) {
          trace_accum[d] += trace_local[d];
        }
      }
    }
  }

  if (rc_error != 0) {
    Rcpp::stop(
      "backend api->trace_hinv_t failed for shard %d (code %d: %s)",
      rc_group,
      rc_error,
      trace_hinv_t_strerror(rc_error)
    );
  }

  Rcpp::NumericVector trace_out(Nparams);
  for (size_t d = 0; d < Nparams; ++d) {
    trace_out[static_cast<R_xlen_t>(d)] = trace_accum[d];
  }

  return trace_out;
}
