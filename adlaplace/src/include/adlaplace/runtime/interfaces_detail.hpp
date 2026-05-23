#ifndef ADLAPLACE_INTERFACES_DETAIL_HPP
#define ADLAPLACE_INTERFACES_DETAIL_HPP

// adlaplace internal: resolve handles, sparsity shards, shard index lists.

#include <Rcpp.h>
#include <vector>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/runtime/ad_fun_pack.hpp"

adlaplace_adpack_handle* shard_handle(ad_fun* backend, size_t shard);

SEXP ad_fun_handle_sexp(const Rcpp::List& ad_fun_list);

ad_fun* ad_fun_from_list(const Rcpp::List& ad_fun_list);

ad_fun* ad_fun_from_handle(SEXP handle);

ad_fun* resolve_ad_fun_eval(SEXP ad_fun_ptr);

ad_fun* resolve_ad_fun_laplace(const Rcpp::S4& ad_fun_s4);

std::vector<size_t> resolve_shard_indices(
  size_t n_shards,
  const Rcpp::IntegerVector& shards);

Rcpp::List sparsity_shard_from_handle(adlaplace_adpack_handle* h);

#endif
