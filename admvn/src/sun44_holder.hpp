#ifndef ADMVN_SUN44_HOLDER_HPP
#define ADMVN_SUN44_HOLDER_HPP

#include "sun44_tape.hpp"

#include <Rcpp.h>
#include <vector>

namespace admvn {

class Sun44TapeHolder {
public:
  Sun44TapeBundle bundle;
  std::vector<double> seed_par;

  Sun44TapeHolder(Sun44TapeBundle b, std::vector<double> par)
    : bundle(std::move(b)), seed_par(std::move(par)) {}
};

inline void sun44_holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<Sun44TapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

inline Sun44TapeHolder* sun44_holder_from_sexp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn SUN44 tape pointer");
  }
  return static_cast<Sun44TapeHolder*>(R_ExternalPtrAddr(ptr));
}

}  // namespace admvn

#endif
