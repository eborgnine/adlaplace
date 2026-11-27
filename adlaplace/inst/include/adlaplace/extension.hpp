#ifndef ADLAPLACE_EXTENSION_HPP
#define ADLAPLACE_EXTENSION_HPP

// Public API for backend packages and single-file sourceCpp models that
// compile custom log densities and CppAD tape construction in their own shared
// library.
//
// On macOS, tape recording and evaluation must live in the same .so as the
// density functions. Backends do not link adlaplace.so; registration helpers
// are header-only and must be called from the one translation unit that
// includes register_impl.hpp (where make_ad_pack_ptr and packs_to_ad_fun
// are defined). Rcpp attributes are parsed before macro expansion, so the
// export and ADLAPLACE_DEFINE_BACKEND stay in the user's file.
//
// Single-file sourceCpp model. sourcecpp.hpp includes this header plus the
// once-only eval and registration implementations. RcppEigen is required
// (backend.hpp includes <Eigen/SparseCore>); the cpp17 plugin matters when
// R is older than 4.3:
//
//   // [[Rcpp::depends(adlaplace, RCppAD, RcppEigen)]]
//   // [[Rcpp::plugins(cpp17)]]
//   #include "adlaplace/sourcecpp.hpp"
//
//   CppAD::vector<CppAD::AD<double>> my_density(
//     const CppAD::vector<CppAD::AD<double>>& par,
//     const density_data& data, const Config& config) { ... }
//
//   ADLAPLACE_DEFINE_BACKEND(my_make_shard)
//
//   // [[Rcpp::export]]
//   SEXP my_ad_pack(SEXP data, Rcpp::List config) {
//     return parameters_ad_pack_ptr(data, config, my_density, my_make_shard);
//   }
//
// observations_ad_pack_ptr is the same shape for an observation density
// (the function takes size_t Dgroup; one tape per obs_groups column).
// A backend package instead exports get_ad_pack_raw_obs /
// get_ad_pack_raw_parameters; random densities export create_ad_shard_<name>
// looked up via data@package.

#include <cstddef>
#include <vector>

#include <Rcpp.h>

#include "adlaplace/backend.hpp"
#include "adlaplace/density_data.hpp"
#include "adlaplace/ad_pack.hpp"
#include "adlaplace/eta.hpp"
#include "adlaplace/rviews.hpp"
#include "adlaplace/eval.hpp"

#if defined(__GNUC__)
#define ADLAPLACE_EXTENSION_EXPORT __attribute__((visibility("default")))
#else
#define ADLAPLACE_EXTENSION_EXPORT
#endif

SEXP make_ad_pack_ptr(ad_pack* groups);

ad_pack* packs_to_ad_fun(
  std::vector<AdTape>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  ShardFactory factory);

// Tape a parameters density (one tape over all mapped beta/gamma and all
// theta) and return an ad_pack_ptr. Pass config through unchanged.
inline SEXP parameters_ad_pack_ptr(
    SEXP data, const Rcpp::List& config,
    LogDensSingleDataFn log_dens, ShardFactory factory) {
  const density_data model(data);
  std::vector<AdTape> tapes;
  tapes.push_back(build_ad_fun_parameters(model, config, log_dens));
  return make_ad_pack_ptr(packs_to_ad_fun(
      std::move(tapes), model.num_beta, model.num_theta, factory));
}

// Same for an observation density: one tape per obs_groups column.
// log_dens takes the group index. Pass config through unchanged
// (config$num_threads selects the parallel build).
inline SEXP observations_ad_pack_ptr(
    SEXP data, const Rcpp::List& config,
    LogDensObsFn log_dens, ShardFactory factory) {
  const density_data model(data);
  std::vector<AdTape> tapes = build_ad_fun_obs(model, config, log_dens);
  return make_ad_pack_ptr(packs_to_ad_fun(
      std::move(tapes), model.num_beta, model.num_theta, factory));
}

// Include eval_impl.hpp exactly once in the backend .cpp before this macro.
#define ADLAPLACE_DEFINE_BACKEND(factory_sym) \
  ADLAPLACE_EXTENSION_EXPORT ad_shard* factory_sym(AdTape&& pack) { \
    return new EvalShard(std::move(pack), factory_sym); \
  }

#endif
