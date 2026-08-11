#ifndef ADMVN_SUN43_HOLDER_HPP
#define ADMVN_SUN43_HOLDER_HPP

#include "sun43_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class Sun43TapeHolder {
public:
  Sun43TapeBundle bundle;
  std::vector<double> seed_par;

  Sun43TapeHolder(Sun43TapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun43_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<Sun43TapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline Sun43TapeHolder* sun43_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN43 tape pointer");
  }
  return static_cast<Sun43TapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
