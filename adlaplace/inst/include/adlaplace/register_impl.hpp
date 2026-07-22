// Include in exactly one .cpp per shared library (adlaplace.so or a backend .so),
// alongside eval_impl.hpp. Provides packs_to_ad_fun / make_ad_pack_ptr so backends
// do not need to link adlaplace.so.

#include "adlaplace/extension.hpp"

#include <R.h>
#include <Rinternals.h>

void ad_fun_destroy(ad_pack* groups) {
  if (!groups) return;
  for (ad_shard* shard : groups->fun) {
    delete shard;
  }
  delete groups;
}

void adfun_finalizer(SEXP ext) {
  ad_pack* groups = static_cast<ad_pack*>(R_ExternalPtrAddr(ext));
  ad_fun_destroy(groups);
  R_ClearExternalPtr(ext);
}

SEXP make_ad_pack_ptr(ad_pack* groups) {
  SEXP handle = R_MakeExternalPtr(static_cast<void*>(groups), R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(handle, adfun_finalizer, TRUE);
  Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("ad_pack_ptr"));
  return handle;
}

ad_pack* packs_to_ad_fun(
  std::vector<AdTape>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  ShardFactory factory) {

  if (!factory) {
    Rcpp::stop("packs_to_ad_fun: factory is NULL");
  }

  auto* groups = new ad_pack();
  groups->abi_version = ADLAPLACE_ABI_VERSION;
  groups->fun.reserve(packs.size());
  for (size_t g = 0; g < packs.size(); ++g) {
    packs[g].shard_index = g;
    packs[g].n_beta = n_beta;
    packs[g].n_theta = n_theta;
    groups->fun.push_back(factory(std::move(packs[g])));
  }
  return groups;
}
