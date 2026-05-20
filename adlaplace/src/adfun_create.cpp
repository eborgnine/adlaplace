#include <Rcpp.h>

#include <Rinternals.h>
#include "adlaplace/runtime/interfaces.hpp"

//' Build AD pack external pointer handle
//'
//' Constructs grouped CppAD tapes and returns an opaque \code{ad_groups}
//' handle (fun only; attach Hessian templates via \code{adlaplace_attach_hessian()}).
//'
//' @param data Model data list passed to the backend builder.
//' @param config Model configuration list passed to the backend builder.
//'
//' @return External pointer of class \code{adlaplace_handle_ptr}.
//'
//' @seealso \code{\link{getAdFun}}
//' @export
// [[Rcpp::export]]
SEXP build_adfun(Rcpp::List data, Rcpp::List config) {
  return build_adfun_h(data, config);
}

//' Attach hessian_map() result to an ad_groups handle
//'
//' Copies outer/inner templates and maps into \code{ad_groups}.
//'
//' @param handle External pointer from \code{build_adfun()} or \code{getAdFun()}.
//' @param hessian_pack List returned by \code{hessian_map()}.
//' @export
// [[Rcpp::export]]
void adlaplace_attach_hessian(SEXP handle, Rcpp::List hessian_map) {
  ad_groups* groups = ad_groups_from_handle(handle);
  ad_groups_attach_hessians_from_list(*groups, hessian_map);
  ad_groups_attach_chol_pattern(*groups, hessian_map);
}

//' @title C++ backend entry points
//' @name adlaplace_cpp
//' @description Low-level C++ entry points exposed to R via Rcpp.
//'
//' @section Sign convention:
//' \code{jointLogDens()}, \code{grad()}, and \code{hessian()} evaluate the
//' \strong{joint log density} \eqn{\ell(x)} and its derivatives (maximization sign).
//' \code{all_derivs()} and the objective inside \code{inner_opt()} use the
//' \strong{negative} log density \eqn{-\ell(x)} and its derivatives for
//' \pkg{trustOptim} minimization. At the same \code{x} and \code{ad_fun},
//' \code{all_derivs()$fval == -jointLogDens(ad_fun, x)},
//' \code{all_derivs()$gradient == -grad(ad_fun, x)}, and
//' \code{all_derivs()$hessian == -hessian(ad_fun, x)} (outer, full parameter vector).
//' @param handle External pointer returned by \code{build_adfun()}.
//' @param group 0-based group index for per-shard sparsity queries.
//' @param data Model data list required by the backend builder.
//' @param config Model configuration list required by the backend builder.
//' @param x Numeric parameter vector of length \code{Nparams}.
//' @param inner Logical scalar for inner-\eqn{\gamma} vs outer derivatives.
//' @param Sgroups Optional integer vector of 0-based group indices.
//' @param LinvPt,LinvPtColumns,verbose,num_threads See \code{traceHinvT()}.
//'
//' @rdname adlaplace_cpp
// [[Rcpp::export]]
int adlaplace_n_groups(SEXP handle) {
  ad_groups* groups = ad_groups_from_handle(handle);
  return static_cast<int>(groups->fun.size());
}

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
Rcpp::List adlaplace_get_sparse_sizes(SEXP handle, int group) {
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

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
Rcpp::List adlaplace_get_sparse_pattern(SEXP handle, int group) {
  ad_groups* groups = ad_groups_from_handle(handle);
  return sparsity_shard_from_handle(shard_handle(groups, static_cast<size_t>(group)));
}
