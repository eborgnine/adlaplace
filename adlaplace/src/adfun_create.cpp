#include <Rcpp.h>

#include <Rinternals.h>
#include "adlaplace/api/register.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"

//' Build raw AD handle for observation shards only
//'
//' @param model An \code{ad_model} S4 object.
//' @param config Model configuration list.
//' @param name Registered observation density name (e.g. \code{"neg_binom_obs"}).
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_obs(SEXP model, Rcpp::List config, std::string name) {
  ad_fun* groups = get_ad_fun_raw_obs_h(model, config, name);
  return make_ad_fun_ptr(groups);
}

//' Build raw AD handle for a random-effect shard
//'
//' @param model An \code{ad_model} S4 object (maps for the term).
//' @param precision Precision list (\code{Q} vector for \code{random_diagonal}).
//' @param config Model configuration list.
//' @param name Registered random density name (e.g. \code{"random_diagonal"}).
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_random(SEXP model, Rcpp::List precision, Rcpp::List config, std::string name) {
  ad_fun* groups = get_ad_fun_raw_random_h(model, precision, config, name);
  return make_ad_fun_ptr(groups);
}

//' Build raw AD handle for a parameters shard
//'
//' @param model An \code{ad_model} S4 object.
//' @param config Model configuration list.
//' @param name Registered parameters density name (e.g. \code{"neg_binom_extra"}).
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_parameters(SEXP model, Rcpp::List config, std::string name) {
  ad_fun* groups = get_ad_fun_raw_parameters_h(model, config, name);
  return make_ad_fun_ptr(groups);
}

//' Merge partial AD handles into one raw handle
//'
//' Concatenates shards from \code{\link{ad_fun_ptr}} (or other partial
//' builders) in list order. Does not attach \code{hessian_map}; use
//' \code{\link{ad_fun}} when templates are needed.
//'
//' @param handles List of external pointers (\code{ad_fun_ptr}).
//' @return Combined external pointer of class \code{ad_fun_ptr}.
//' @seealso \code{\link{ad_fun_ptr}}, \code{\link{ad_fun}}
//' @keywords internal
// [[Rcpp::export(name = c_ad_fun_ptr)]]
SEXP c_ad_fun_ptr(Rcpp::List handles) {
  std::vector<ad_fun*> parts;
  parts.reserve(handles.size());
  for (int i = 0; i < handles.size(); ++i) {
    SEXP h = handles[i];
    if (TYPEOF(h) != EXTPTRSXP) {
      Rcpp::stop("handles[[%d]] must be an external pointer", i + 1);
    }
    ad_fun* g = static_cast<ad_fun*>(R_ExternalPtrAddr(h));
    if (!g) {
      Rcpp::stop("handles[[%d]] external pointer is NULL", i + 1);
    }
    parts.push_back(g);
    R_ClearExternalPtr(h);
  }
  ad_fun* merged = combine_ad_fun(parts);
  return make_ad_fun_ptr(merged);
}

//' Attach hessian_map() result to an ad_fun handle
//'
//' Copies outer/inner templates and maps into the C++ \code{ad_fun} handle.
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @param hessian_pack List returned by \code{hessian_map()}.
//' @export
// [[Rcpp::export]]
void adlaplace_attach_hessian(SEXP handle, Rcpp::List hessian_pack) {
  ad_fun* groups = ad_fun_from_handle(handle);
  ad_fun_attach_hessians_from_list(*groups, hessian_pack);
}

//' Number of AD shards in an \code{ad_fun_ptr} handle
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @return Integer count of groups (shards).
//' @export
// [[Rcpp::export]]
int n_groups(SEXP handle) {
  ad_fun* groups = ad_fun_from_handle(handle);
  return static_cast<int>(groups->fun.size());
}

//' Sparse structure sizes for one AD shard
//'
//' Layout and sparsity sizes for one AD shard
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @param group 0-based group index.
//' @return List with \code{n_inner}, \code{n_outer}, \code{n_beta}, \code{n_theta},
//'   and \code{nnz_grad_*}, \code{nnz_hes_*}.
//' @export
// [[Rcpp::export]]
Rcpp::List get_sizes(SEXP handle, int group) {
  ad_fun* groups = ad_fun_from_handle(handle);
  adlaplace_adpack_handle* h = shard_handle(groups, static_cast<size_t>(group));
  if (!h->api->get_sizes) {
    Rcpp::stop("backend api->get_sizes is NULL");
  }

  int n_inner = 0;
  int n_outer = 0;
  int n_beta = 0;
  int n_theta = 0;
  int nnz_grad_inner = 0;
  int nnz_grad_outer = 0;
  int nnz_hes_inner = 0;
  int nnz_hes_outer = 0;
  if (h->api->get_sizes(
      h->ctx, &n_inner, &n_outer, &n_beta, &n_theta,
      &nnz_grad_inner, &nnz_grad_outer,
      &nnz_hes_inner, &nnz_hes_outer) != 0) {
    Rcpp::stop("backend api->get_sizes failed for group %d", group);
  }

  return Rcpp::List::create(
    Rcpp::Named("n_inner") = n_inner,
    Rcpp::Named("n_outer") = n_outer,
    Rcpp::Named("n_beta") = n_beta,
    Rcpp::Named("n_theta") = n_theta,
    Rcpp::Named("nnz_grad_inner") = nnz_grad_inner,
    Rcpp::Named("nnz_grad_outer") = nnz_grad_outer,
    Rcpp::Named("nnz_hes_inner") = nnz_hes_inner,
    Rcpp::Named("nnz_hes_outer") = nnz_hes_outer
  );
}

//' Sparse index patterns for one AD shard
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @param group 0-based group index.
//' @return List with \code{grad}, \code{grad_inner}, \code{row_hess}, \code{col_hess}, etc.
//' @export
// [[Rcpp::export]]
Rcpp::List get_sparse_pattern(SEXP handle, int group) {
  ad_fun* groups = ad_fun_from_handle(handle);
  return sparsity_shard_from_handle(shard_handle(groups, static_cast<size_t>(group)));
}
