#ifndef ADLAPLACE_REGISTER_HPP
#define ADLAPLACE_REGISTER_HPP

// Backend registration: export raw ad_groups handle to R.
// Include from Rinterfaces.cpp after model tapes are defined in objectiveFunction.cpp.

#include <Rcpp.h>
#include <memory>
#include <vector>

#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/api/backend.hpp"

// Defined in adlaplace (src/adpack_api.cpp); link backend package to adlaplace.
extern "C" const adlaplace_adpack_api adlaplace_AD_API;

std::vector<GroupPack> get_ad_fun(const Data& data, const Config& config);

inline void backend_destroy(void* vctx) {
  delete static_cast<GroupPack*>(vctx);
}

inline void ad_groups_destroy(ad_groups* groups) {
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

inline void adgroups_finalizer(SEXP ext) {
  ad_groups* groups = static_cast<ad_groups*>(R_ExternalPtrAddr(ext));
  ad_groups_destroy(groups);
  R_ClearExternalPtr(ext);
}

inline ad_groups* get_ad_fun_raw_h(
  const Data& data,
  const Config& config) {

  std::vector<GroupPack> packs = get_ad_fun(data, config);
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

inline SEXP get_ad_fun_raw_h(
  const Rcpp::List& data,
  const Rcpp::List& config) {

  const Data dataC(data);
  const Config configC(config);

  ad_groups* groups = get_ad_fun_raw_h(dataC, configC);

  SEXP handle = R_MakeExternalPtr(static_cast<void*>(groups), R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(handle, adgroups_finalizer, TRUE);
  Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));
  return handle;
}

#endif
