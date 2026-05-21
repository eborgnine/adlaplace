#ifndef ADLAPLACE_INTERFACES_DETAIL_HPP
#define ADLAPLACE_INTERFACES_DETAIL_HPP

// adlaplace internal: resolve handles, sparsity shards, group lists.

#include <Rcpp.h>
#include <vector>

#include "adlaplace/api/register.hpp"
#include "adlaplace/runtime/ad_groups_pack.hpp"

Rcpp::List extract_sparsity(const std::vector<GroupPack>& ad_fun);

static inline adlaplace_adpack_handle* shard_handle(ad_groups* groups, size_t g) {
  if (!groups) Rcpp::stop("ad_groups handle is NULL");
  if (g >= groups->fun.size()) {
    Rcpp::stop("group index %d out of range [0, %d]", (int)g, (int)groups->fun.size() - 1);
  }
  adlaplace_adpack_handle* h = groups->fun[g];
  if (!h) Rcpp::stop("ad_groups.fun[%d] is NULL", (int)g);
  if (!h->api) Rcpp::stop("ad_groups.fun[%d] has NULL api", (int)g);
  if (!h->ctx) Rcpp::stop("ad_groups.fun[%d] has NULL ctx", (int)g);
  if (!h->api->f) Rcpp::stop("ad_groups.fun[%d] api->f is NULL", (int)g);
  return h;
}

static inline SEXP ad_fun_handle_sexp(const Rcpp::List& ad_fun) {
  if (ad_fun.containsElementNamed("ad_fun")) {
    return ad_fun["ad_fun"];
  }
  Rcpp::stop("ad_fun list must contain component 'ad_fun'");
}

static inline ad_groups* get_ad_groups(const Rcpp::List& ad_fun) {
  SEXP handle_ptr = ad_fun_handle_sexp(ad_fun);
  auto* groups = static_cast<ad_groups*>(R_ExternalPtrAddr(handle_ptr));
  if (!groups) {
    Rcpp::stop("ad_groups external pointer is NULL (cleared?)");
  }

  ad_groups_attach_hessians_from_list(*groups, ad_fun);
  ad_groups_attach_chol_pattern(*groups, ad_fun);
  return groups;
}

static inline ad_groups* ad_groups_from_handle(SEXP handle) {
  auto* groups = static_cast<ad_groups*>(R_ExternalPtrAddr(handle));
  if (!groups) {
    Rcpp::stop("ad_groups external pointer is NULL (cleared?)");
  }
  return groups;
}

static inline ad_groups* resolve_ad_groups(SEXP ad_fun) {
  if (Rf_isNull(ad_fun)) {
    Rcpp::stop("ad_fun must not be NULL");
  }
  if (TYPEOF(ad_fun) == EXTPTRSXP) {
    return ad_groups_from_handle(ad_fun);
  }
  if (TYPEOF(ad_fun) == VECSXP) {
    return get_ad_groups(Rcpp::as<Rcpp::List>(ad_fun));
  }
  Rcpp::stop("ad_fun must be a list from get_ad_fun() or an external pointer from get_ad_fun_raw()");
}

static inline std::vector<size_t> resolve_groups(
  size_t Ngroups,
  const Rcpp::IntegerVector& Sgroups) {
  std::vector<size_t> groups;
  if (Sgroups.size() == 0) {
    groups.resize(Ngroups);
    for (size_t g = 0; g < Ngroups; ++g) groups[g] = g;
    return groups;
  }

  groups.reserve(Sgroups.size());
  for (R_xlen_t k = 0; k < Sgroups.size(); ++k) {
    if (Sgroups[k] == NA_INTEGER) {
      Rcpp::stop("Sgroups contains NA at position %d", (int)k + 1);
    }
    if (Sgroups[k] < 0 || static_cast<size_t>(Sgroups[k]) >= Ngroups) {
      Rcpp::stop("Sgroups index %d out of range [0, %d]", Sgroups[k], (int)Ngroups - 1);
    }
    groups.push_back(static_cast<size_t>(Sgroups[k]));
  }
  return groups;
}

inline Rcpp::List sparsity_shard_from_handle(adlaplace_adpack_handle* h) {
  if (!h->api->get_sparse_sizes) {
    Rcpp::stop("backend api->get_sparse_sizes is NULL");
  }
  if (!h->api->get_sparse_pattern) {
    Rcpp::stop("backend api->get_sparse_pattern is NULL");
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
    Rcpp::stop("backend api->get_sparse_sizes failed");
  }

  std::vector<int> grad_inner(static_cast<size_t>(nnz_grad_inner));
  std::vector<int> grad_outer(static_cast<size_t>(nnz_grad_outer));
  std::vector<int> row_hes_inner(static_cast<size_t>(nnz_hes_inner));
  std::vector<int> col_hes_inner(static_cast<size_t>(nnz_hes_inner));
  std::vector<int> row_hes_outer(static_cast<size_t>(nnz_hes_outer));
  std::vector<int> col_hes_outer(static_cast<size_t>(nnz_hes_outer));

  if (h->api->get_sparse_pattern(
        h->ctx,
        grad_inner.data(), grad_outer.data(),
        row_hes_inner.data(), col_hes_inner.data(),
        row_hes_outer.data(), col_hes_outer.data()) != 0) {
    Rcpp::stop("backend api->get_sparse_pattern failed");
  }

  return Rcpp::List::create(
    Rcpp::Named("grad") = Rcpp::IntegerVector(
      grad_outer.begin(), grad_outer.end()),
    Rcpp::Named("grad_inner") = Rcpp::IntegerVector(
      grad_inner.begin(), grad_inner.end()),
    Rcpp::Named("row_hess") = Rcpp::IntegerVector(
      row_hes_outer.begin(), row_hes_outer.end()),
    Rcpp::Named("col_hess") = Rcpp::IntegerVector(
      col_hes_outer.begin(), col_hes_outer.end()),
    Rcpp::Named("row_hess_inner") = Rcpp::IntegerVector(
      row_hes_inner.begin(), row_hes_inner.end()),
    Rcpp::Named("col_hess_inner") = Rcpp::IntegerVector(
      col_hes_inner.begin(), col_hes_inner.end())
  );
}

#endif
