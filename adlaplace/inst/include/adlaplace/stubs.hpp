#ifndef ADLAPLACE_STUBS_HPP
#define ADLAPLACE_STUBS_HPP

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include "adlaplace/backend.hpp"

typedef void (*adlaplace_hessian_pack_from_list_c_t)(SEXP, HessianPack*, HessianPack*);

inline std::array<HessianPack, 2> adlaplace_hessianPackFromList(const Rcpp::List& x) {
  static adlaplace_hessian_pack_from_list_c_t fn = nullptr;
  if (fn == nullptr) {
    fn = reinterpret_cast<adlaplace_hessian_pack_from_list_c_t>(
      R_GetCCallable("adlaplace", "adlaplace_hessian_pack_from_list_c")
    );
    if (fn == nullptr) {
      Rcpp::stop("Unable to resolve callable adlaplace_hessian_pack_from_list_c");
    }
  }

  std::array<HessianPack, 2> out;
  fn(x, &out[0], &out[1]);
  return out;
}

#endif
