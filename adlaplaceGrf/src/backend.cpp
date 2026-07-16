// Backend boilerplate for adlaplaceGrf (FEM Matérn random densities).
// Taping and evaluation must live in this DSO (macOS CppAD requirement).

#include <Rcpp.h>

#include "adlaplace/extension.hpp"

#include "adlaplace/eval_impl.hpp"

ADLAPLACE_DEFINE_BACKEND(adlaplace_grf_make_shard)
