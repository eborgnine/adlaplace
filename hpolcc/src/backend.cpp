#include <Rcpp.h>

#include "adlaplace/adfun.hpp"
#include "adlaplace/extension.hpp"

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config,
  const size_t Dgroup);

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Config& config);

#include "adlaplace/eval_impl.hpp"
#include "adlaplace/register_impl.hpp"

ADLAPLACE_DEFINE_BACKEND(hpolcc_make_shard)

static LogDensObsFn resolve_obs_density(const std::string& name) {
  if (name == "dirichlet_multinomial") return logDensObs;
  Rcpp::stop("unknown observation density: %s", name.c_str());
}

static LogDensSingleDataFn resolve_extra_density(const std::string& name) {
  if (name == "dirichlet_multinomial_extra") return logDensExtra;
  Rcpp::stop("unknown parameters density: %s", name.c_str());
}

//' Build observation-shard \code{ad_fun_ptr} in this package's shared library.
//'
//' @param model \code{ad_data} S4 object.
//' @param config Model configuration list.
//' @param name Observation density name (e.g. \code{"dirichlet_multinomial"}).
//' @return External pointer of class \code{ad_fun_ptr}.
// [[Rcpp::export]]
SEXP get_ad_fun_raw_obs(SEXP model, Rcpp::List config, std::string name) {
  const ad_data ad_model(model);
  std::vector<GroupPack> packs = build_ad_fun_obs(
    ad_model, config, resolve_obs_density(name));
  ad_fun* groups = packs_to_ad_fun(
    std::move(packs),
    ad_model.num_beta,
    ad_model.num_theta,
    hpolcc_make_shard);
  return make_ad_fun_ptr(groups);
}

//' Build parameters-shard \code{ad_fun_ptr} in this package's shared library.
//'
//' @param model \code{ad_data} S4 object.
//' @param config Model configuration list.
//' @param name Parameters density name (e.g. \code{"dirichlet_multinomial_extra"}).
//' @return External pointer of class \code{ad_fun_ptr}.
// [[Rcpp::export]]
SEXP get_ad_fun_raw_parameters(SEXP model, Rcpp::List config, std::string name) {
  const ad_data ad_model(model);
  GroupPack pack = build_ad_fun_parameters(
    ad_model, config, resolve_extra_density(name));
  std::vector<GroupPack> packs;
  packs.push_back(std::move(pack));
  ad_fun* groups = packs_to_ad_fun(
    std::move(packs),
    ad_model.num_beta,
    ad_model.num_theta,
    hpolcc_make_shard);
  return make_ad_fun_ptr(groups);
}
