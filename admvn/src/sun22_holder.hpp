#ifndef ADMVN_SUN22_HOLDER_HPP
#define ADMVN_SUN22_HOLDER_HPP

#include "sun22_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class Sun22TapeHolder {
public:
  Sun22TapeBundle bundle;
  std::vector<double> seed_par;

  Sun22TapeHolder(Sun22TapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun22_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<Sun22TapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline Sun22TapeHolder* sun22_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN22 tape pointer");
  }
  return static_cast<Sun22TapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
