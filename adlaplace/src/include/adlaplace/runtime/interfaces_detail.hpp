#ifndef ADLAPLACE_INTERFACES_DETAIL_HPP
#define ADLAPLACE_INTERFACES_DETAIL_HPP

// adlaplace internal: resolve handles, sparsity shards, group lists.

#include <Rcpp.h>
#include <vector>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/runtime/ad_groups_pack.hpp"

adlaplace_adpack_handle* shard_handle(ad_groups* groups, size_t g);

SEXP ad_fun_handle_sexp(const Rcpp::List& ad_fun);

ad_groups* get_ad_groups(const Rcpp::List& ad_fun);

ad_groups* ad_groups_from_handle(SEXP handle);

ad_groups* resolve_ad_groups(SEXP ad_fun);

std::vector<size_t> resolve_groups(
  size_t Ngroups,
  const Rcpp::IntegerVector& Sgroups);

Rcpp::List sparsity_shard_from_handle(adlaplace_adpack_handle* h);

#endif
