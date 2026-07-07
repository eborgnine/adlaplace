#include "adlaplace/api/register.hpp"
#include "adlaplace/extension/adfun_pack.hpp"

#include "adlaplace/creators/adfun_obs.hpp"
#include "adlaplace/creators/adfun_single.hpp"
#include "adlaplace/densities/binomial.hpp"
#include "adlaplace/densities/gaussian.hpp"
#include "adlaplace/densities/nbinom.hpp"
#include "adlaplace/densities/poisson.hpp"
#include "adlaplace/densities/random_diagonal.hpp"
#include "adlaplace/densities/random_mult.hpp"

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

static LogDensSingleRandomDiagFn resolve_random_density(const std::string& name) {
  if (name == "random_diagonal" || name == "random") return random_diagonal;
  if (name == "random_mult") return random_mult;
  Rcpp::stop("unknown random density: %s", name.c_str());
}

void backend_destroy(void* vctx) {
  delete static_cast<GroupPack*>(vctx);
}

void ad_fun_destroy(ad_fun* groups) {
  if (!groups) return;
  for (adlaplace_adpack_handle* h : groups->fun) {
    if (!h) continue;
    if (h->api && h->api->destroy && h->ctx) {
      h->api->destroy(h->ctx);
    }
    delete h;
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

static ad_fun* packs_to_ad_fun_layout(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  const adlaplace_adpack_api* api) {

  if (!api) {
    Rcpp::stop("adlaplace_packs_to_ad_fun: api is NULL");
  }

  auto* groups = new ad_fun();
  groups->fun.reserve(packs.size());
  for (size_t g = 0; g < packs.size(); ++g) {
    auto* pack = new GroupPack(std::move(packs[g]));
    pack->shard_index = g;
    pack->n_beta = n_beta;
    pack->n_theta = n_theta;
    auto* h = new adlaplace_adpack_handle();
    h->api = api;
    h->ctx = static_cast<void*>(pack);
    groups->fun.push_back(h);
  }
  return groups;
}

ad_fun* packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta) {
  return packs_to_ad_fun_layout(std::move(packs), n_beta, n_theta, &adlaplace_AD_API);
}

ad_fun* packs_to_ad_fun_api(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  const adlaplace_adpack_api* api) {
  return packs_to_ad_fun_layout(std::move(packs), n_beta, n_theta, api);
}

SEXP adlaplace_make_ad_fun_ptr(ad_fun* groups) {
  return make_ad_fun_ptr(groups);
}

ad_fun* adlaplace_packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta) {
  return packs_to_ad_fun(std::move(packs), n_beta, n_theta);
}

ad_fun* adlaplace_packs_to_ad_fun_api(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  const adlaplace_adpack_api* api) {
  return packs_to_ad_fun_api(std::move(packs), n_beta, n_theta, api);
}

static GroupPack* first_group_pack(ad_fun* groups) {
  if (!groups || groups->fun.empty()) {
    Rcpp::stop("ad_fun handle has no shards");
  }
  for (adlaplace_adpack_handle* h : groups->fun) {
    if (h && h->ctx) {
      return static_cast<GroupPack*>(h->ctx);
    }
  }
  Rcpp::stop("ad_fun handle has no valid shards");
}

ad_fun* combine_ad_fun(const std::vector<ad_fun*>& parts) {
  if (parts.empty()) {
    Rcpp::stop("c_ad_fun_ptr requires at least one handle");
  }
  GroupPack* layout = first_group_pack(parts[0]);
  const std::size_t n_beta = layout->n_beta;
  const std::size_t n_theta = layout->n_theta;

  std::vector<GroupPack> merged;
  std::vector<const adlaplace_adpack_api*> merged_apis;
  size_t shard_index = 0;

  for (ad_fun* part : parts) {
    if (!part) continue;
    for (adlaplace_adpack_handle* h : part->fun) {
      if (!h || !h->ctx) continue;
      GroupPack* pack = static_cast<GroupPack*>(h->ctx);
      pack->shard_index = shard_index++;
      pack->n_beta = n_beta;
      pack->n_theta = n_theta;
      merged_apis.push_back(h->api ? h->api : &adlaplace_AD_API);
      merged.push_back(std::move(*pack));
      delete pack;
      h->ctx = nullptr;
    }
    part->fun.clear();
    ad_fun_destroy(part);
  }

  auto* groups = new ad_fun();
  groups->fun.reserve(merged.size());
  for (size_t g = 0; g < merged.size(); ++g) {
    auto* pack = new GroupPack(std::move(merged[g]));
    auto* h = new adlaplace_adpack_handle();
    h->api = merged_apis[g];
    h->ctx = static_cast<void*>(pack);
    groups->fun.push_back(h);
  }
  return groups;
}

ad_fun* get_ad_fun_raw_obs_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& obs_name) {

  const ad_data model(model_sexp);
  std::vector<GroupPack> packs = build_ad_fun_obs(model, config, resolve_obs_density(obs_name));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta);
}

ad_fun* get_ad_fun_raw_random_h(
  SEXP model_sexp,
  SEXP precision,
  const Rcpp::List& config,
  const std::string& single_name) {

  const ad_data model(model_sexp);
  GroupPack pack = build_ad_fun_random(model, precision, config, resolve_random_density(single_name));
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta);
}

ad_fun* get_ad_fun_raw_parameters_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& single_name) {

  const ad_data model(model_sexp);
  GroupPack pack = build_ad_fun_parameters(model, config, resolve_extra_density(single_name));
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta);
}
