#include <Rcpp.h>

#include <Rinternals.h>
#include "adlaplace/register.hpp"
#include "adlaplace/runtime.hpp"

//' Build raw AD handle for observation shards only
//'
//' @param model An \code{ad_data} S4 object.
//' @param config Model configuration list.
//' @param name Registered observation density name (e.g. \code{"nbinom_obs"}).
//' @return External pointer of class \code{ad_fun_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_obs(SEXP model, Rcpp::List config, std::string name) {
  ad_fun* groups = get_ad_fun_raw_obs_h(model, config, name);
  return make_ad_fun_ptr(groups);
}

//' Build raw AD handle for a parameters shard
//'
//' @param model An \code{ad_data} S4 object.
//' @param config Model configuration list.
//' @param name Registered parameters density name (e.g. \code{"nbinom_extra"}).
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
//' @keywords internal
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
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List get_sizes(SEXP handle, int group) {
  ad_fun* groups = ad_fun_from_handle(handle);
  adlaplace_shard* shard = shard_handle(groups, static_cast<size_t>(group));

  int n_inner = 0;
  int n_outer = 0;
  int n_beta = 0;
  int n_theta = 0;
  int nnz_grad_inner = 0;
  int nnz_grad_outer = 0;
  int nnz_hes_inner = 0;
  int nnz_hes_outer = 0;
  if (shard->get_sizes(
      &n_inner, &n_outer, &n_beta, &n_theta,
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
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List get_sparse_pattern(SEXP handle, int group) {
  ad_fun* groups = ad_fun_from_handle(handle);
  return sparsity_shard_from_handle(shard_handle(groups, static_cast<size_t>(group)));
}

//' Build-time owner thread for one AD shard
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @param group 0-based group index.
//' @return Integer owner thread id recorded when the shard was prewarmed.
//' @keywords internal
// [[Rcpp::export]]
int get_thread_owner(SEXP handle, int group) {
  ad_fun* groups = ad_fun_from_handle(handle);
  adlaplace_shard* shard = shard_handle(groups, static_cast<size_t>(group));
  return static_cast<int>(shard->pack.owner_thread);
}

//' Configured \code{num_threads} from \code{ad_fun()} thread assignment
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @return Integer \code{num_threads} if \code{ad_fun()} assigned threads; otherwise \code{NA}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::IntegerVector get_configured_num_threads(SEXP handle) {
  ad_fun* groups = ad_fun_from_handle(handle);
  if (!groups->num_threads_configured) {
    return Rcpp::IntegerVector(NA_INTEGER);
  }
  return Rcpp::IntegerVector::create(
    static_cast<int>(groups->configured_num_threads)
  );
}

//' Whether a shard has an assigned OpenMP owner thread
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @param group 0-based group index.
//' @return Logical scalar.
//' @keywords internal
// [[Rcpp::export]]
bool get_owner_thread_assigned(SEXP handle, int group) {
  ad_fun* groups = ad_fun_from_handle(handle);
  adlaplace_shard* shard = shard_handle(groups, static_cast<size_t>(group));
  return shard->pack.owner_thread_assigned;
}

//' Assign OpenMP owner threads to all shards
//'
//' Sets \code{owner_thread = shard_index \% num_threads} on every shard.
//' Called from \code{ad_fun()} after \code{c()}. Required before parallel
//' \code{inner_opt} / \code{trace_hinv_t}.
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @param num_threads Positive integer thread count.
//' @keywords internal
// [[Rcpp::export]]
void assign_owner_threads(SEXP handle, int num_threads) {
  if (num_threads < 1) {
    Rcpp::stop("num_threads must be a positive integer");
  }
  ad_fun* groups = ad_fun_from_handle(handle);
  const std::size_t n_threads = static_cast<std::size_t>(num_threads);
  for (std::size_t s = 0; s < groups->fun.size(); ++s) {
    adlaplace_shard* shard = shard_handle(groups, s);
    shard->pack.owner_thread = shard->pack.shard_index % n_threads;
    shard->pack.owner_thread_assigned = true;
  }
  groups->configured_num_threads = n_threads;
  groups->num_threads_configured = true;
}

//' Deep copy of an \code{ad_fun_ptr} handle
//'
//' Clones CppAD tapes and sparsity patterns into a new external pointer.
//' The source handle is unchanged (unlike \code{c()} on \code{ad_fun_ptr}, which moves
//' shards and clears sources).
//'
//' @param handle External pointer of class \code{ad_fun_ptr}.
//' @return New \code{ad_fun_ptr} with independent C++ state.
//' @keywords internal
// [[Rcpp::export(name = "clone_ad_fun_ptr_")]]
SEXP clone_ad_fun_ptr_impl(SEXP handle) {
  ad_fun* src = ad_fun_from_handle(handle);
  ad_fun* copy = clone_ad_fun(src);
  return make_ad_fun_ptr(copy);
}
