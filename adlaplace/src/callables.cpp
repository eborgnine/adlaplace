#include <Rcpp.h>
#include <R_ext/Rdynload.h>

#include "adlaplace/rviews.hpp"

extern "C" SEXP adlaplace_grad_c(SEXP xSEXP, SEXP backendContextSEXP, SEXP innerSEXP, SEXP SgroupsSEXP);
extern "C" SEXP adlaplace_hess_c(SEXP xSEXP, SEXP backendContextSEXP, SEXP innerSEXP, SEXP SgroupsSEXP);
extern "C" SEXP adlaplace_jointLogDens_c(SEXP xSEXP, SEXP backendContextSEXP, SEXP SgroupsSEXP);
extern "C" SEXP adlaplace_traceHinvT_c(
  SEXP xSEXP,
  SEXP backendContextSEXP,
  SEXP LinvPtSEXP,
  SEXP LinvPtColumnsSEXP,
  SEXP SgroupsSEXP
);

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
    R_RegisterCCallable(
      "adlaplace",
      "adlaplace_grad_c",
      reinterpret_cast<DL_FUNC>(adlaplace_grad_c)
    );
    R_RegisterCCallable(
      "adlaplace",
      "adlaplace_hess_c",
      reinterpret_cast<DL_FUNC>(adlaplace_hess_c)
    );
    R_RegisterCCallable(
      "adlaplace",
      "adlaplace_jointLogDens_c",
      reinterpret_cast<DL_FUNC>(adlaplace_jointLogDens_c)
    );
    R_RegisterCCallable(
      "adlaplace",
      "adlaplace_traceHinvT_c",
      reinterpret_cast<DL_FUNC>(adlaplace_traceHinvT_c)
    );
    registered = true;
  }
  return R_NilValue;
}
