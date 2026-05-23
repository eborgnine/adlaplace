#ifndef ADLAPLACE_REGISTER_HPP
#define ADLAPLACE_REGISTER_HPP

// adlaplace internal: build and merge ad_fun handles for R export.

#include <Rcpp.h>
#include <memory>
#include <string>
#include <vector>

#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/creators/ad_model.hpp"
#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/api/backend.hpp"

extern "C" const adlaplace_adpack_api adlaplace_AD_API;

void register_adlaplace_default_densities();

void backend_destroy(void* vctx);
void ad_fun_destroy(ad_fun* groups);
void adfun_finalizer(SEXP ext);
SEXP make_ad_fun_ptr(ad_fun* groups);

ad_fun* packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta);

ad_fun* combine_ad_fun(const std::vector<ad_fun*>& parts);

ad_fun* get_ad_fun_raw_obs_h(
  SEXP model,
  const Rcpp::List& config,
  const std::string& obs_name);

ad_fun* get_ad_fun_raw_random_h(
  SEXP model,
  const Rcpp::List& precision,
  const Rcpp::List& config,
  const std::string& single_name);

ad_fun* get_ad_fun_raw_parameters_h(
  SEXP model,
  const Rcpp::List& config,
  const std::string& single_name);

#endif
