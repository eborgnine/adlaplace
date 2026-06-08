#include <Rcpp.h>

#include "adlaplace/api/log_dens_fn.hpp"
#include "adlaplace/creators/adfun_obs.hpp"
#include "adlaplace/creators/adfun_single.hpp"
#include "adlaplace/extension/adfun_pack.hpp"
#include "adlaplace/extension/backend_adpack_api.hpp"

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup);

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list);

void backend_destroy(void* vctx) {
  delete static_cast<GroupPack*>(vctx);
}

#include "adlaplace/extension/backend_eval_amalgam.hpp"

ADLAPLACE_DEFINE_BACKEND_ADPACK_API(adlaplace_example_AD_API)

static LogDensObsFn resolve_obs_density(const std::string& name) {
  if (name == "skewnormal_obs") return logDensObs;
  Rcpp::stop("unknown observation density: %s", name.c_str());
}

static LogDensSingleDataFn resolve_extra_density(const std::string& name) {
  if (name == "skewnormal_extra") return logDensExtra;
  Rcpp::stop("unknown parameters density: %s", name.c_str());
}

//' Build observation-shard \code{ad_fun_ptr} in this package's shared library.
//'
//' Records CppAD tapes in \code{adlaplaceExample.so} (same DSO as custom densities).
//' @param model \code{ad_data} S4 object.
//' @param config Model configuration list.
//' @param name Observation density name (e.g. \code{"skewnormal_obs"}).
//' @return External pointer of class \code{ad_fun_ptr}.
// [[Rcpp::export]]
SEXP get_ad_fun_raw_obs(SEXP model, Rcpp::List config, std::string name) {
  const ad_data ad_model(model);
  std::vector<GroupPack> packs = build_ad_fun_obs(
    ad_model, config, resolve_obs_density(name));
  ad_fun* groups = adlaplace_packs_to_ad_fun_api(
    std::move(packs),
    ad_model.num_beta,
    ad_model.num_theta,
    &adlaplace_example_AD_API);
  return adlaplace_make_ad_fun_ptr(groups);
}

//' Build parameters-shard \code{ad_fun_ptr} in this package's shared library.
//'
//' @param model \code{ad_data} S4 object.
//' @param config Model configuration list.
//' @param name Parameters density name (e.g. \code{"skewnormal_extra"}).
//' @return External pointer of class \code{ad_fun_ptr}.
// [[Rcpp::export]]
SEXP get_ad_fun_raw_parameters(SEXP model, Rcpp::List config, std::string name) {
  const ad_data ad_model(model);
  GroupPack pack = build_ad_fun_parameters(
    ad_model, config, resolve_extra_density(name));
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  ad_fun* groups = adlaplace_packs_to_ad_fun_api(
    std::move(packs),
    ad_model.num_beta,
    ad_model.num_theta,
    &adlaplace_example_AD_API);
  return adlaplace_make_ad_fun_ptr(groups);
}
