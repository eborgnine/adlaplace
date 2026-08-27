#include "adlaplace/register.hpp"
#include "adlaplace/extension.hpp"
#include "adlaplace/ad_pack.hpp"
#include "adlaplace/densities.hpp"
#include "adlaplace/runtime.hpp"

#define ADLAPLACE_MATH_LGAMMA_DEFINE
#define ADLAPLACE_MATH_LOG_ERFC_DEFINE
#include "adlaplace/atomics.hpp"

#include "adlaplace/register_impl.hpp"

extern ad_shard* adlaplace_make_ad_shard(AdTape&&);

static LogDensObsFn resolve_obs_density(const std::string& name) {
  if (name == "nbinom_obs") return nbinom_obs;
  if (name == "poisson_obs") return poisson_obs;
  if (name == "gaussian_obs") return gaussian_obs;
  if (name == "binomial_obs") return binomial_obs;
  if (name == "dirichlet_multinomial") return dirichlet_multinomial;
  Rcpp::stop("unknown observation density: %s", name.c_str());
}

static LogDensSingleDataFn resolve_extra_density(const std::string& name) {
  if (name == "nbinom_extra") return nbinom_extra;
  if (name == "gaussian_extra") return gaussian_extra;
  if (name == "dirichlet_multinomial_extra") return dirichlet_multinomial_extra;
  if (name == "exp_prior") return exp_prior;
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

static ad_shard* first_shard(ad_pack* groups) {
  if (!groups || groups->fun.empty()) {
    Rcpp::stop("ad_pack handle has no shards");
  }
  for (ad_shard* shard : groups->fun) {
    if (shard) {
      return shard;
    }
  }
  Rcpp::stop("ad_pack handle has no valid shards");
}

ad_pack* combine_ad_fun(const std::vector<ad_pack*>& parts) {
  if (parts.empty()) {
    Rcpp::stop("c_ad_pack_ptr requires at least one handle");
  }
  ad_shard* layout = first_shard(parts[0]);
  const std::size_t n_beta = layout->pack.n_beta;
  const std::size_t n_theta = layout->pack.n_theta;

  std::vector<AdTape> merged;
  std::vector<ShardFactory> merged_factories;
  size_t shard_index = 0;

  for (ad_pack* part : parts) {
    if (!part) continue;
    if (part->abi_version != ADLAPLACE_ABI_VERSION) {
      Rcpp::stop(
        "ad_pack ABI mismatch (got %d, need %d): rebuild the extension package against this adlaplace",
        part->abi_version,
        ADLAPLACE_ABI_VERSION
      );
    }
    for (ad_shard* shard : part->fun) {
      if (!shard) continue;
      shard->pack.shard_index = shard_index++;
      shard->pack.n_beta = n_beta;
      shard->pack.n_theta = n_theta;
      merged_factories.push_back(shard->factory ? shard->factory : adlaplace_make_ad_shard);
      merged.push_back(std::move(shard->pack));
      delete shard;
    }
    part->fun.clear();
    ad_fun_destroy(part);
  }

  auto* groups = new ad_pack();
  groups->abi_version = ADLAPLACE_ABI_VERSION;
  groups->fun.reserve(merged.size());
  for (size_t g = 0; g < merged.size(); ++g) {
    groups->fun.push_back(merged_factories[g](std::move(merged[g])));
  }
  return groups;
}

ad_pack* get_ad_pack_raw_obs_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& obs_name) {

  const density_data model(model_sexp);
  std::vector<AdTape> packs = build_ad_fun_obs(model, config, resolve_obs_density(obs_name));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta, adlaplace_make_ad_shard);
}

ad_pack* get_ad_pack_raw_parameters_h(
  SEXP model_sexp,
  const Rcpp::List& config,
  const std::string& single_name) {

  const density_data model(model_sexp);
  AdTape pack = build_ad_fun_parameters(model, config, resolve_extra_density(single_name));
  std::vector<AdTape> packs;
  packs.push_back(std::move(pack));
  return packs_to_ad_fun(std::move(packs), model.num_beta, model.num_theta, adlaplace_make_ad_shard);
}
