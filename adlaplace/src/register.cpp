#include "adlaplace/api/register.hpp"

#include "adlaplace/creators/adfun_obs.hpp"
#include "adlaplace/creators/adfun_single.hpp"

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
  std::size_t n_theta) {

  auto* groups = new ad_fun();
  groups->fun.reserve(packs.size());
  for (size_t g = 0; g < packs.size(); ++g) {
    auto* pack = new GroupPack(std::move(packs[g]));
    pack->shard_index = g;
    pack->n_beta = n_beta;
    pack->n_theta = n_theta;
    auto* h = new adlaplace_adpack_handle();
    h->api = &adlaplace_AD_API;
    h->ctx = static_cast<void*>(pack);
    groups->fun.push_back(h);
  }
  return groups;
}

ad_fun* packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta) {
  return packs_to_ad_fun_layout(std::move(packs), n_beta, n_theta);
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
  size_t shard_index = 0;

  for (ad_fun* part : parts) {
    if (!part) continue;
    for (adlaplace_adpack_handle* h : part->fun) {
      if (!h || !h->ctx) continue;
      GroupPack* pack = static_cast<GroupPack*>(h->ctx);
      pack->shard_index = shard_index++;
      pack->n_beta = n_beta;
      pack->n_theta = n_theta;
      merged.push_back(std::move(*pack));
      delete pack;
      h->ctx = nullptr;
    }
    part->fun.clear();
    ad_fun_destroy(part);
  }

  return packs_to_ad_fun_layout(std::move(merged), n_beta, n_theta);
}

ad_fun* get_ad_fun_raw_obs_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& obs_name) {

  register_adlaplace_default_densities();
  const ad_model model(model_sexp);
  std::vector<GroupPack> packs = build_ad_fun_obs(model, config, obs_name);
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta);
}

ad_fun* get_ad_fun_raw_random_h(
  SEXP model_sexp,
  const Rcpp::List& precision,
  const Rcpp::List& config,
  const std::string& single_name) {

  register_adlaplace_default_densities();
  const ad_model model(model_sexp);
  GroupPack pack = build_ad_fun_random(model, precision, config, single_name);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta);
}

ad_fun* get_ad_fun_raw_parameters_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& single_name) {

  register_adlaplace_default_densities();
  const ad_model model(model_sexp);
  GroupPack pack = build_ad_fun_parameters(model, config, single_name);
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta);
}
