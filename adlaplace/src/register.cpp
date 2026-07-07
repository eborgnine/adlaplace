#include "adlaplace/register.hpp"
#include "adlaplace/extension.hpp"
#include "adlaplace/adfun.hpp"
#include "adlaplace/densities.hpp"
#include "adlaplace/runtime.hpp"

#define ADLAPLACE_MATH_LGAMMA_DEFINE
#define ADLAPLACE_MATH_LOG_ERFC_DEFINE
#include "adlaplace/atomics.hpp"

extern adlaplace_shard* adlaplace_make_shard(GroupPack&&);

static LogDensObsFn resolve_obs_density(const std::string& name) {
  if (name == "nbinom_obs") return nbinom_obs;
  if (name == "poisson_obs") return poisson_obs;
  if (name == "gaussian_obs") return gaussian_obs;
  if (name == "binomial_obs") return binomial_obs;
  Rcpp::stop("unknown observation density: %s", name.c_str());
}

static LogDensSingleDataFn resolve_extra_density(const std::string& name) {
  if (name == "nbinom_extra") return nbinom_extra;
  if (name == "gaussian_extra") return gaussian_extra;
  Rcpp::stop("unknown parameters density: %s", name.c_str());
}

void adlaplace_init_atomics() {
  static bool initialized = false;
  if (initialized) return;
  initialized = true;

  adlaplace_debug_print_load_banner();

  adlaplace_init_lgamma_atomic();
  adlaplace_init_log_erfc_atomic();
}

void ad_fun_destroy(ad_fun* groups) {
  if (!groups) return;
  for (adlaplace_shard* shard : groups->fun) {
    delete shard;
  }
  delete groups;
}

void adfun_finalizer(SEXP ext) {
  ad_fun* groups = static_cast<ad_fun*>(R_ExternalPtrAddr(ext));
  ad_fun_destroy(groups);
  R_ClearExternalPtr(ext);
}

SEXP make_ad_fun_ptr(ad_fun* groups) {
  SEXP handle = R_MakeExternalPtr(static_cast<void*>(groups), R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(handle, adfun_finalizer, TRUE);
  Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("ad_fun_ptr"));
  return handle;
}

ad_fun* packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  ShardFactory factory) {

  if (!factory) {
    Rcpp::stop("packs_to_ad_fun: factory is NULL");
  }

  auto* groups = new ad_fun();
  groups->fun.reserve(packs.size());
  for (size_t g = 0; g < packs.size(); ++g) {
    packs[g].shard_index = g;
    packs[g].n_beta = n_beta;
    packs[g].n_theta = n_theta;
    groups->fun.push_back(factory(std::move(packs[g])));
  }
  return groups;
}

static adlaplace_shard* first_shard(ad_fun* groups) {
  if (!groups || groups->fun.empty()) {
    Rcpp::stop("ad_fun handle has no shards");
  }
  for (adlaplace_shard* shard : groups->fun) {
    if (shard) {
      return shard;
    }
  }
  Rcpp::stop("ad_fun handle has no valid shards");
}

ad_fun* combine_ad_fun(const std::vector<ad_fun*>& parts) {
  if (parts.empty()) {
    Rcpp::stop("c_ad_fun_ptr requires at least one handle");
  }
  adlaplace_shard* layout = first_shard(parts[0]);
  const std::size_t n_beta = layout->pack.n_beta;
  const std::size_t n_theta = layout->pack.n_theta;

  std::vector<GroupPack> merged;
  std::vector<ShardFactory> merged_factories;
  size_t shard_index = 0;

  for (ad_fun* part : parts) {
    if (!part) continue;
    for (adlaplace_shard* shard : part->fun) {
      if (!shard) continue;
      shard->pack.shard_index = shard_index++;
      shard->pack.n_beta = n_beta;
      shard->pack.n_theta = n_theta;
      merged_factories.push_back(shard->factory ? shard->factory : adlaplace_make_shard);
      merged.push_back(std::move(shard->pack));
      delete shard;
    }
    part->fun.clear();
    ad_fun_destroy(part);
  }

  auto* groups = new ad_fun();
  groups->fun.reserve(merged.size());
  for (size_t g = 0; g < merged.size(); ++g) {
    groups->fun.push_back(merged_factories[g](std::move(merged[g])));
  }
  return groups;
}

ad_fun* get_ad_fun_raw_obs_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& obs_name) {

  const ad_data model(model_sexp);
  std::vector<GroupPack> packs = build_ad_fun_obs(model, config, resolve_obs_density(obs_name));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta, adlaplace_make_shard);
}

ad_fun* get_ad_fun_raw_parameters_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& single_name) {

  const ad_data model(model_sexp);
  GroupPack pack = build_ad_fun_parameters(model, config, resolve_extra_density(single_name));
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta, adlaplace_make_shard);
}
