#ifndef ADLAPLACE_EXTENSION_ADFUN_PACK_HPP
#define ADLAPLACE_EXTENSION_ADFUN_PACK_HPP

// Public API for backend packages (e.g. adlaplaceExample) that compile custom
// log densities and CppAD tape construction in their own shared library.
//
// On macOS, build_ad_fun_obs / build_ad_fun_parameters must be compiled into
// the same .so as the density functions; registering pointers alone is not enough.

#include <cstddef>
#include <vector>

#include <Rcpp.h>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/api/backend.hpp"
#include "adlaplace/api/log_dens_fn.hpp"
#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/creators/adfun_obs.hpp"
#include "adlaplace/creators/adfun_single.hpp"
#include "adlaplace/creators/rviews.hpp"

#if defined(__GNUC__)
#define ADLAPLACE_EXTENSION_EXPORT __attribute__((visibility("default")))
#else
#define ADLAPLACE_EXTENSION_EXPORT
#endif

ADLAPLACE_EXTENSION_EXPORT SEXP adlaplace_make_ad_fun_ptr(ad_fun* groups);

ADLAPLACE_EXTENSION_EXPORT ad_fun* adlaplace_packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta);

ADLAPLACE_EXTENSION_EXPORT ad_fun* adlaplace_packs_to_ad_fun_api(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  const adlaplace_adpack_api* api);

#endif
