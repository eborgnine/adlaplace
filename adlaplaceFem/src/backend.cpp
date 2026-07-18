// Backend boilerplate for adlaplaceFem (FEM Matérn random densities).
// Taping and evaluation must live in this DSO (macOS CppAD requirement).
// register_impl.hpp provides packs_to_ad_fun / make_ad_fun_ptr (no adlaplace.so link).

#include <Rcpp.h>

#include "adlaplace/extension.hpp"

#include "adlaplace/eval_impl.hpp"
#include "adlaplace/register_impl.hpp"

ADLAPLACE_DEFINE_BACKEND(adlaplace_fem_make_shard)
