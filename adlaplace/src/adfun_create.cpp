#include <Rcpp.h>

#include <Rinternals.h>
#include "adlaplace/runtime/interfaces.hpp"

//' Build raw AD pack external pointer handle
//'
//' Constructs grouped CppAD tapes and returns an opaque \code{ad_groups}
//' handle (fun only; attach Hessian templates via \code{adlaplace_attach_hessian()}).
//'
//' @param data Model data list passed to the backend builder.
//' @param config Model configuration list passed to the backend builder.
//'
//' @return External pointer of class \code{adlaplace_handle_ptr}.
//'
//' @seealso \code{\link{get_ad_fun}}
//' @export
// [[Rcpp::export]]
SEXP get_ad_fun_raw(Rcpp::List data, Rcpp::List config) {
  return get_ad_fun_raw_h(data, config);
}

//' Attach hessian_map() result to an ad_groups handle
//'
//' Copies outer/inner templates and maps into \code{ad_groups}.
//'
//' @param handle External pointer from \code{get_ad_fun_raw()} or \code{get_ad_fun()}.
//' @param hessian_pack List returned by \code{hessian_map()}.
//' @export
// [[Rcpp::export]]
void adlaplace_attach_hessian(SEXP handle, Rcpp::List hessian_map) {
  ad_groups* groups = ad_groups_from_handle(handle);
  ad_groups_attach_hessians_from_list(*groups, hessian_map);
  ad_groups_attach_chol_pattern(*groups, hessian_map);
}

//' Number of AD shards in an \code{ad_groups} handle
//'
//' @param handle External pointer from \code{get_ad_fun_raw()} or \code{get_ad_fun()}.
//' @return Integer count of groups (shards).
//' @export
// [[Rcpp::export]]
int n_groups(SEXP handle) {
  ad_groups* groups = ad_groups_from_handle(handle);
  return static_cast<int>(groups->fun.size());
}

//' Sparse structure sizes for one AD shard
//'
//' @param handle External pointer from \code{get_ad_fun_raw()} or \code{get_ad_fun()}.
//' @param group 0-based group index.
//' @return List with \code{n_inner}, \code{n_outer}, \code{nnz_grad_*}, \code{nnz_hes_*}.
//' @export
// [[Rcpp::export]]
Rcpp::List get_sparse_sizes(SEXP handle, int group) {
  ad_groups* groups = ad_groups_from_handle(handle);
  adlaplace_adpack_handle* h = shard_handle(groups, static_cast<size_t>(group));
  if (!h->api->get_sparse_sizes) {
    Rcpp::stop("backend api->get_sparse_sizes is NULL");
  }

  int n_inner = 0;
  int n_outer = 0;
  int nnz_grad_inner = 0;
  int nnz_grad_outer = 0;
  int nnz_hes_inner = 0;
  int nnz_hes_outer = 0;
  if (h->api->get_sparse_sizes(
      h->ctx, &n_inner, &n_outer,
      &nnz_grad_inner, &nnz_grad_outer,
      &nnz_hes_inner, &nnz_hes_outer) != 0) {
    Rcpp::stop("backend api->get_sparse_sizes failed for group %d", group);
  }

  return Rcpp::List::create(
    Rcpp::Named("n_inner") = n_inner,
    Rcpp::Named("n_outer") = n_outer,
    Rcpp::Named("nnz_grad_inner") = nnz_grad_inner,
    Rcpp::Named("nnz_grad_outer") = nnz_grad_outer,
    Rcpp::Named("nnz_hes_inner") = nnz_hes_inner,
    Rcpp::Named("nnz_hes_outer") = nnz_hes_outer
  );
}

//' Sparse index patterns for one AD shard
//'
//' @param handle External pointer from \code{get_ad_fun_raw()} or \code{get_ad_fun()}.
//' @param group 0-based group index.
//' @return List with \code{grad}, \code{grad_inner}, \code{row_hess}, \code{col_hess}, etc.
//' @export
// [[Rcpp::export]]
Rcpp::List get_sparse_pattern(SEXP handle, int group) {
  ad_groups* groups = ad_groups_from_handle(handle);
  return sparsity_shard_from_handle(shard_handle(groups, static_cast<size_t>(group)));
}
