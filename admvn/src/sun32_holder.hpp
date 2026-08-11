#ifndef ADMVN_SUN32_HOLDER_HPP
#define ADMVN_SUN32_HOLDER_HPP

#include "sun32_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class Sun32TapeHolder {
public:
  Sun32TapeBundle bundle;
  std::vector<double> seed_par;

  Sun32TapeHolder(Sun32TapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun32_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<Sun32TapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline Sun32TapeHolder* sun32_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN32 tape pointer");
  }
  return static_cast<Sun32TapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
