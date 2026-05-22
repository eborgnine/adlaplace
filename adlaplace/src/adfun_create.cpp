#include <Rcpp.h>

#include <Rinternals.h>
#include "adlaplace/api/register.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"
#include "chol_update.hpp"

//' Build raw AD handle for observation shards only
//'
//' @param data Model data list.
//' @param config Model configuration list.
//' @param name Registered observation density name (e.g. \code{"neg_binom_obs"}).
//' @return External pointer of class \code{adlaplace_handle_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_obs(Rcpp::List data, Rcpp::List config, std::string name) {
  ad_groups* groups = get_ad_fun_raw_obs_h(data, config, name);
  return make_ad_groups_handle(groups);
}

//' Build raw AD handle for one single-density shard
//'
//' @param data Model data list.
//' @param config Model configuration list.
//' @param name Registered single density name (e.g. \code{"random_diagonal"}).
//' @return External pointer of class \code{adlaplace_handle_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_single(Rcpp::List data, Rcpp::List config, std::string name) {
  ad_groups* groups = get_ad_fun_raw_single_h(data, config, name);
  return make_ad_groups_handle(groups);
}

//' Merge partial AD handles into one raw handle
//'
//' Concatenates shards from \code{\link{get_ad_fun_raw}} (or other partial
//' builders) in list order. Does not attach \code{hessian_map}; use
//' \code{\link{get_ad_fun}} when templates are needed.
//'
//' @param handles List of external pointers (\code{adlaplace_handle_ptr}).
//' @param config Model configuration list (sets \code{n_beta}/\code{n_gamma} on shards).
//' @return Combined external pointer of class \code{adlaplace_handle_ptr}.
//' @seealso \code{\link{get_ad_fun_raw}}, \code{\link{get_ad_fun}}
//' @export
// [[Rcpp::export]]
SEXP combine(Rcpp::List handles, Rcpp::List config) {
  const Config configC(config);
  std::vector<ad_groups*> parts;
  parts.reserve(handles.size());
  for (int i = 0; i < handles.size(); ++i) {
    SEXP h = handles[i];
    if (TYPEOF(h) != EXTPTRSXP) {
      Rcpp::stop("handles[[%d]] must be an external pointer", i + 1);
    }
    ad_groups* g = static_cast<ad_groups*>(R_ExternalPtrAddr(h));
    if (!g) {
      Rcpp::stop("handles[[%d]] external pointer is NULL", i + 1);
    }
    parts.push_back(g);
    R_ClearExternalPtr(h);
  }
  ad_groups* merged = combine_ad_groups(parts, configC);
  return make_ad_groups_handle(merged);
}

//' Attach hessian_map() result to an ad_groups handle
//'
//' Copies outer/inner templates and maps into \code{ad_groups}.
//'
//' @param handle External pointer from \code{get_ad_fun_raw()} or \code{get_ad_fun()}.
//' @param hessian_pack List returned by \code{hessian_map()}.
//' @export
// [[Rcpp::export]]
void adlaplace_attach_hessian(SEXP handle, Rcpp::List hessian_pack) {
  ad_groups* groups = ad_groups_from_handle(handle);
  ad_groups_attach_hessians_from_list(*groups, hessian_pack);
  ad_groups_attach_chol_pattern(*groups, hessian_pack);
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
