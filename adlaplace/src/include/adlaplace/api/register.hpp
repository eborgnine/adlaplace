#ifndef ADLAPLACE_REGISTER_HPP
#define ADLAPLACE_REGISTER_HPP

// adlaplace internal: build and merge ad_groups handles for R export.

#include <Rcpp.h>
#include <memory>
#include <string>
#include <vector>

#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/api/backend.hpp"

extern "C" const adlaplace_adpack_api adlaplace_AD_API;

void register_adlaplace_default_densities();

void backend_destroy(void* vctx);
void ad_groups_destroy(ad_groups* groups);
void adgroups_finalizer(SEXP ext);
SEXP make_ad_groups_handle(ad_groups* groups);

ad_groups* packs_to_ad_groups(
  std::vector<GroupPack>&& packs,
  const Config& config);

ad_groups* combine_ad_groups(
  const std::vector<ad_groups*>& parts,
  const Config& config);

ad_groups* get_ad_fun_raw_obs_h(
  const Rcpp::List& data,
  const Rcpp::List& config,
  const std::string& obs_name);

ad_groups* get_ad_fun_raw_single_h(
  const Rcpp::List& data,
  const Rcpp::List& config,
  const std::string& single_name);

#endif
