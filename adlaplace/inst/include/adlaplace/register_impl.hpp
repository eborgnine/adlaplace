// Include in exactly one .cpp per shared library (adlaplace.so or a backend .so),
// alongside eval_impl.hpp. Provides packs_to_ad_fun / make_ad_fun_ptr so backends
// do not need to link adlaplace.so.

#include "adlaplace/extension.hpp"

#include <R.h>
#include <Rinternals.h>

void ad_fun_destroy(ad_fun* groups) {
  if (!groups) return;
  for (adlaplace_shard* shard : groups->fun) {
    delete shard;
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

ad_fun* packs_to_ad_fun(
  std::vector<GroupPack>&& packs,
  std::size_t n_beta,
  std::size_t n_theta,
  ShardFactory factory) {

  if (!factory) {
    Rcpp::stop("packs_to_ad_fun: factory is NULL");
  }

  auto* groups = new ad_fun();
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
