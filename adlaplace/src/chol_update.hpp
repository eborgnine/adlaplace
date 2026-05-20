#ifndef ADLAPLACE_CHOL_UPDATE_HPP
#define ADLAPLACE_CHOL_UPDATE_HPP

#include <Rcpp.h>

#include "adlaplace/runtime/backend.hpp"

CholPattern chol_pattern_from_sexp(SEXP chol_inner);

#endif
