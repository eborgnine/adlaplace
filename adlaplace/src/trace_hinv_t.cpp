#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <vector>

#include "adlaplace/ompad.hpp"
#include "adlaplace/runtime/interfaces.hpp"

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  const int num_threads,
  SEXP Sgroups
) {
  adlaplace_adpack_handle* h = get_handle(ad_fun);
  if (!h->api->trace_hinv_t) {
    Rcpp::stop("ad_fun api->trace_hinv_t is NULL");
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
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
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  int rc_error = 0;
  int rc_group = -1;
  const int n_groups = static_cast<int>(groups.size());

  if (n_groups > 0) {
    cppad_parallel_setup(static_cast<std::size_t>(num_threads));

#pragma omp parallel num_threads(num_threads)
    {
      std::vector<double> trace_local(Nparams, 0.0);
      int rc_local = 0;
      int rc_group_local = -1;

#pragma omp for schedule(static,1)
      for (int gidx = 0; gidx < n_groups; ++gidx) {
        const int gi = static_cast<int>(groups[static_cast<std::size_t>(gidx)]);
        const int rc = h->api->trace_hinv_t(
          h->ctx,
          &gi,
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
          rc_group_local = gi;
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
    Rcpp::stop("backend api->trace_hinv_t failed for group %d with code %d", rc_group, rc_error);
  }

  Rcpp::NumericVector trace_out(Nparams);
  for (size_t d = 0; d < Nparams; ++d) {
    trace_out[static_cast<R_xlen_t>(d)] = trace_accum[d];
  }

  return trace_out;
}
