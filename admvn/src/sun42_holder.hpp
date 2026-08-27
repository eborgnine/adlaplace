#ifndef ADMVN_SUN42_HOLDER_HPP
#define ADMVN_SUN42_HOLDER_HPP

#include "sun42_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class Sun42TapeHolder {
public:
  Sun42TapeBundle bundle;
  std::vector<double> seed_par;

  Sun42TapeHolder(Sun42TapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun42_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<Sun42TapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline Sun42TapeHolder* sun42_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN42 tape pointer");
  }
  return static_cast<Sun42TapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
