#ifndef ADLAPLACE_REGISTER_HPP
#define ADLAPLACE_REGISTER_HPP

#include <Rcpp.h>
#include <string>
#include <vector>

#include "adlaplace/rviews.hpp"
#include "adlaplace/ad_data.hpp"
#include "adlaplace/backend.hpp"
#include "adlaplace/extension.hpp"

void adlaplace_init_atomics();

void ad_fun_destroy(ad_fun* groups);
void adfun_finalizer(SEXP ext);
SEXP make_ad_fun_ptr(ad_fun* groups);

ad_fun* packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  ShardFactory factory);

ad_fun* combine_ad_fun(const std::vector<ad_fun*>& parts);

ad_fun* get_ad_fun_raw_obs_h(
  SEXP model,
  const Rcpp::List& config,
  const std::string& obs_name);

ad_fun* get_ad_fun_raw_parameters_h(
  SEXP model,
  const Rcpp::List& config,
  const std::string& single_name);

#endif
