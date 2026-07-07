#include <Rcpp.h>

#include <vector>

#include "adlaplace/runtime.hpp"
#include "adlaplace/trace_hinv_t_runtime.hpp"

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

  const std::vector<double> trace_accum = adlaplace_trace::trace_hinv_t_impl(
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
