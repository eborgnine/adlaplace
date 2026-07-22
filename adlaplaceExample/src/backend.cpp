// backend.cpp -- the boilerplate every adlaplace backend package needs.
//
// A backend package supplies model-specific log densities (see
// objectiveFunction.cpp) and records their CppAD tapes inside its own
// shared library; adlaplace supplies the Laplace approximation, inner
// optimization, derivatives, and parallel evaluation. The required pieces,
// in order, are:
//
//   1. declarations of your log-density functions,
//   2. #include "adlaplace/eval_impl.hpp" and register_impl.hpp -- exactly once, so CppAD tape
//      *evaluation* runs in this DSO (macOS requires taping and evaluation
//      in the same shared object),
//   3. ADLAPLACE_DEFINE_BACKEND(<unique_factory_name>) -- exports a shard
//      factory adlaplace calls through virtual dispatch,
//   4. resolve_*_density() maps -- from the ad_pack name on the R side
//      (e.g. "skewnormal_obs") to your C++ function pointers,
//   5. Rcpp exports get_ad_pack_raw_obs() / get_ad_pack_raw_parameters();
//      adlaplace::ad_pack_ptr() dispatches to these via the `package` slot
//      on the model term. Random-effect densities use create_ad_shard_<name>
//      in adlaplace (or another package) instead — no backend random export.
//
// To write a new backend: copy this file, rename the density names and the
// factory symbol, and point the resolvers at your own functions.

#include <Rcpp.h>

#include "adlaplace/ad_pack.hpp"
#include "adlaplace/extension.hpp"

CppAD::vector<CppAD::AD<double>> logDensObs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const density_data& model,
  const Config& config,
  const size_t Dgroup);

CppAD::vector<CppAD::AD<double>> logDensExtra(
  const CppAD::vector<CppAD::AD<double>>& x,
  const density_data& model,
  const Config& config);

#include "adlaplace/eval_impl.hpp"
#include "adlaplace/register_impl.hpp"

ADLAPLACE_DEFINE_BACKEND(adlaplace_example_make_shard)

static LogDensObsFn resolve_obs_density(const std::string& name) {
  if (name == "skewnormal_obs") return logDensObs;
  Rcpp::stop("unknown observation density: %s", name.c_str());
}

static LogDensSingleDataFn resolve_extra_density(const std::string& name) {
  if (name == "skewnormal_extra") return logDensExtra;
  Rcpp::stop("unknown parameters density: %s", name.c_str());
}

//' Build observation-shard \code{ad_pack_ptr} in this package's shared library.
//'
//' Records CppAD tapes in \code{adlaplaceExample.so} (same DSO as custom densities).
//' @param model \code{density_data} S4 object.
//' @param config Model configuration list.
//' @param name Observation density name (e.g. \code{"skewnormal_obs"}).
//' @return External pointer of class \code{ad_pack_ptr}.
// [[Rcpp::export]]
SEXP get_ad_pack_raw_obs(SEXP model, Rcpp::List config, std::string name) {
  const density_data ad_model(model);
  std::vector<AdTape> packs = build_ad_fun_obs(
    ad_model, config, resolve_obs_density(name));
  ad_pack* groups = packs_to_ad_fun(
    std::move(packs),
    ad_model.num_beta,
    ad_model.num_theta,
    adlaplace_example_make_shard);
  return make_ad_pack_ptr(groups);
}

//' Build parameters-shard \code{ad_pack_ptr} in this package's shared library.
//'
//' @param model \code{density_data} S4 object.
//' @param config Model configuration list.
//' @param name Parameters density name (e.g. \code{"skewnormal_extra"}).
//' @return External pointer of class \code{ad_pack_ptr}.
// [[Rcpp::export]]
SEXP get_ad_pack_raw_parameters(SEXP model, Rcpp::List config, std::string name) {
  const density_data ad_model(model);
  AdTape pack = build_ad_fun_parameters(
    ad_model, config, resolve_extra_density(name));
  std::vector<AdTape> packs;
  packs.push_back(std::move(pack));
  ad_pack* groups = packs_to_ad_fun(
    std::move(packs),
    ad_model.num_beta,
    ad_model.num_theta,
    adlaplace_example_make_shard);
  return make_ad_pack_ptr(groups);
}
