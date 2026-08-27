#ifndef ADMVN_SUN52_HOLDER_HPP
#define ADMVN_SUN52_HOLDER_HPP

#include "sun52_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class Sun52TapeHolder {
public:
  Sun52TapeBundle bundle;
  std::vector<double> seed_par;

  Sun52TapeHolder(Sun52TapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun52_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<Sun52TapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline Sun52TapeHolder* sun52_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN52 tape pointer");
  }
  return static_cast<Sun52TapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
