#include <Rcpp.h>
#include <cppad/cppad.hpp>
#include <omp.h>

#include <vector>

#include "adlaplace/ompad.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"

namespace {

std::vector<std::vector<size_t>> thread_groups_from_backend(ad_fun& backend) {
  const std::size_t n = backend.fun.size();
  std::vector<std::size_t> owners(n);
  std::size_t max_t = 0;
  for (std::size_t s = 0; s < n; ++s) {
    const std::size_t t = pack_ctx(backend.fun[s]->ctx)->owner_thread;
    owners[s] = t;
    if (t > max_t) max_t = t;
  }
  std::vector<std::vector<size_t>> groups(max_t + 1);
  for (std::size_t s = 0; s < n; ++s) {
    groups[owners[s]].push_back(s);
  }
  return groups;
}

std::vector<double> trace_hinv_t_impl(
  ad_fun& backend,
  const std::vector<double>& x,
  const std::vector<int>& LinvPt_p,
  const std::vector<int>& LinvPt_i,
  const std::vector<double>& LinvPt_x,
  const std::size_t LinvPt_ncol,
  const std::vector<int>& LinvPtColumns_p,
  const std::vector<int>& LinvPtColumns_i) {
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

  int rc_error = 0;
  int rc_group = -1;

  const std::vector<std::vector<size_t>> thread_groups = thread_groups_from_backend(backend);
  const int num_threads = static_cast<int>(thread_groups.size());
  cppad_parallel_setup(static_cast<std::size_t>(num_threads));

#pragma omp parallel num_threads(num_threads)
  {
    std::vector<double> trace_local(Nparams, 0.0);
    int rc_local = 0;
    int rc_group_local = -1;
    const std::size_t tid = static_cast<std::size_t>(omp_get_thread_num());
    const std::vector<size_t>& shard_group = thread_groups[tid];

    for (std::size_t s : shard_group) {
      adlaplace_adpack_handle* h = shard_handle(&backend, s);
      if (!h->api->trace_hinv_t) {
        rc_local = 1;
        rc_group_local = static_cast<int>(s);
        continue;
      }
      const int rc = h->api->trace_hinv_t(
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
      for (std::size_t d = 0; d < Nparams; ++d) {
        trace_accum[d] += trace_local[d];
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

  return trace_accum;
}

}  // namespace

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
Rcpp::NumericVector trace_hinv_t(
  SEXP ad_fun_ptr,
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns
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
    LinvPtColumns_i_vec
  );

  return Rcpp::NumericVector(trace_accum.begin(), trace_accum.end());
}
