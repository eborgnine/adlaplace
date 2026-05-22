#include "adlaplace/api/register.hpp"

#include "adlaplace/creators/adfun_obs.hpp"
#include "adlaplace/creators/adfun_single.hpp"

void backend_destroy(void* vctx) {
  delete static_cast<GroupPack*>(vctx);
}

void ad_groups_destroy(ad_groups* groups) {
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

void adgroups_finalizer(SEXP ext) {
  ad_groups* groups = static_cast<ad_groups*>(R_ExternalPtrAddr(ext));
  ad_groups_destroy(groups);
  R_ClearExternalPtr(ext);
}

SEXP make_ad_groups_handle(ad_groups* groups) {
  SEXP handle = R_MakeExternalPtr(static_cast<void*>(groups), R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(handle, adgroups_finalizer, TRUE);
  Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));
  return handle;
}

ad_groups* packs_to_ad_groups(
  std::vector<GroupPack>&& packs,
  const Config& config) {

  auto* groups = new ad_groups();
  groups->fun.reserve(packs.size());
  for (size_t g = 0; g < packs.size(); ++g) {
    auto* pack = new GroupPack(std::move(packs[g]));
    pack->shard_index = g;
    pack->n_beta = config.Nbeta;
    pack->n_gamma = config.Ngamma;
    auto* h = new adlaplace_adpack_handle();
    h->api = &adlaplace_AD_API;
    h->ctx = static_cast<void*>(pack);
    groups->fun.push_back(h);
  }
  return groups;
}

ad_groups* combine_ad_groups(
  const std::vector<ad_groups*>& parts,
  const Config& config) {

  std::vector<GroupPack> merged;
  size_t shard_index = 0;

  for (ad_groups* part : parts) {
    if (!part) continue;
    for (adlaplace_adpack_handle* h : part->fun) {
      if (!h || !h->ctx) continue;
      GroupPack* pack = static_cast<GroupPack*>(h->ctx);
      pack->shard_index = shard_index++;
      pack->n_beta = config.Nbeta;
      pack->n_gamma = config.Ngamma;
      merged.push_back(std::move(*pack));
      delete pack;
      h->ctx = nullptr;
    }
    part->fun.clear();
    ad_groups_destroy(part);
  }

  return packs_to_ad_groups(std::move(merged), config);
}

ad_groups* get_ad_fun_raw_obs_h(
  const Rcpp::List& data,
  const Rcpp::List& config,
  const std::string& obs_name) {

  register_adlaplace_default_densities();
  const Config cfg(config);
  std::vector<GroupPack> packs = build_ad_fun_obs(data, config, obs_name);
  return packs_to_ad_groups(std::move(packs), cfg);
}

ad_groups* get_ad_fun_raw_single_h(
  const Rcpp::List& data,
  const Rcpp::List& config,
  const std::string& single_name) {

  register_adlaplace_default_densities();
  const Config cfg(config);
  GroupPack pack = build_ad_fun_single(data, config, single_name);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_groups(std::move(packs), cfg);
}
