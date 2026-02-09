#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include "adlaplace/rviews.hpp"

extern "C" void adlaplace_hessian_pack_from_list_c(SEXP x, HessianPack* inner, HessianPack* outer) {
  if (inner == nullptr || outer == nullptr) {
    Rcpp::stop("NULL output pointer in adlaplace_hessian_pack_from_list_c");
  }
  std::array<HessianPack, 2> packs = hessianPackFromList(Rcpp::as<Rcpp::List>(x));
  *inner = std::move(packs[0]);
  *outer = std::move(packs[1]);
}

//' Register C-callable entry points
// [[Rcpp::export]]
SEXP register_callables() {
  static bool registered = false;
  if (!registered) {
    R_RegisterCCallable(
      "adlaplace",
      "adlaplace_hessian_pack_from_list_c",
      reinterpret_cast<DL_FUNC>(adlaplace_hessian_pack_from_list_c)
    );
    registered = true;
  }
  return R_NilValue;
}
