#ifndef ADMVN_SUN_HOLDER_HPP
#define ADMVN_SUN_HOLDER_HPP

#include "sun_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class SunTapeHolder {
public:
  SunTapeBundle bundle;
  std::vector<double> seed_par;

  SunTapeHolder(SunTapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<SunTapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline SunTapeHolder* sun_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN tape pointer");
  }
  return static_cast<SunTapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
