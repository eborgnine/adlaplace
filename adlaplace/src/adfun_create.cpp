#include <Rcpp.h>

#include <Rinternals.h>
#include "adlaplace/runtime/interfaces.hpp"

//' C++ backend entry points
//'
//' Low-level C++ entry points exposed to R via Rcpp.
//' These create and operate on backend state returned by \code{getAdFun_r()}.
//' For the default backend this state is a list with an opaque external pointer
//' plus sparsity/Hessian metadata.
//'
//' @param data An R list containing model data objects required by the backend
//'   (used by \code{getAdFun_r()}).
//' @param config An R list of configuration options required by the backend
//'   (used by \code{getAdFun_r()}).
//' @param ad_fun Backend object returned by \code{getAdFun_r()}.
//'   For the default backend this can be either the full returned list
//'   (containing \code{adFun}, \code{sparsity}, \code{hessians}) or the
//'   external pointer in \code{$adFun}.
//' @param x Numeric parameter vector of length \code{Nparams}.
//' @param inner Logical scalar. If \code{TRUE}, evaluate inner-\eqn{\gamma}
//'   derivatives; otherwise evaluate outer/full derivatives.
//' @param Sgroups Optional integer vector of 0-based group indices to evaluate.
//'   If omitted, uses all groups \code{0:(Ngroups-1)}.
//' @param LinvPt Sparse \code{dgCMatrix} for columns of
//'   \eqn{P^\top L^{-1} D^{-1/2}} (or equivalent) used in trace contractions.
//' @param LinvPtColumns Sparse \code{ngCMatrix}/\code{dgCMatrix} mapping
//'   selected columns of \code{LinvPt} to each group.
//' @param verbose Logical scalar indicating whether to print verbose output.
//' @param num_threads Integer specifying the number of threads to use for parallel computation.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun}: external pointer handle to the built AD groups.
//'     Use \code{getAdFun()} in R to also assemble \code{sparsity} and
//'     \code{hessians}.
//'   \item \code{jointLogDens}: scalar objective value summed over groups.
//'   \item \code{grad}: numeric gradient vector.
//'   \item \code{hessian}: symmetric sparse Hessian (\code{dsCMatrix}).
//'   \item \code{traceHinvT}: numeric vector of third-derivative contractions.
//' }
//'
//' @details
//' In the default backend, \code{$adFun} is an opaque external pointer and not
//' user-modifiable. It may hold substantial memory (AD tapes and caches).
//' Do not save backend objects across R sessions.
//'
//' @name adlaplace_cpp

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
SEXP getAdFun(
  Rcpp::List data,
  Rcpp::List config)
{
  return getAdFun_h(data, config);
}

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
SEXP adlaplace_build_groups(Rcpp::List data, Rcpp::List config) {
  return buildAdGroups_h(data, config);
}

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
int adlaplace_n_groups(SEXP handle) {
  adlaplace_adpack_handle* h = get_handle(handle);
  AdGroups* groups = groups_ctx(h->ctx);
  return static_cast<int>(groups->size());
}

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
Rcpp::List adlaplace_get_sizes(SEXP handle) {
  return get_sizes_from_handle(get_handle(handle));
}

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
Rcpp::List adlaplace_get_sparse_sizes(SEXP handle, int group) {
  adlaplace_adpack_handle* h = get_handle(handle);
  if (!h->api->get_sparse_sizes) {
    Rcpp::stop("backend api->get_sparse_sizes is NULL");
  }

  int gi = group;
  int n_inner = 0;
  int n_outer = 0;
  int nnz_grad_inner = 0;
  int nnz_grad_outer = 0;
  int nnz_hes_inner = 0;
  int nnz_hes_outer = 0;
  if (h->api->get_sparse_sizes(
      h->ctx, &gi, &n_inner, &n_outer,
      &nnz_grad_inner, &nnz_grad_outer,
      &nnz_hes_inner, &nnz_hes_outer) != 0) {
    Rcpp::stop("backend api->get_sparse_sizes failed for group %d", gi);
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
  return sparsity_shard_from_handle(get_handle(handle), group);
}

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
void adlaplace_finalize_handle(SEXP handle, Rcpp::List hessians) {
  finalizeAdHandle_h(handle, hessians);
}
