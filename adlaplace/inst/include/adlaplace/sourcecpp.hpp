#ifndef ADLAPLACE_SOURCECPP_HPP
#define ADLAPLACE_SOURCECPP_HPP

// Single-file models (Rcpp::sourceCpp): include once, in the only .cpp.
// eval_impl.hpp and register_impl.hpp define this shared object's tape
// evaluation and ad_pack_ptr builder. A second inclusion duplicates those
// symbols. After this include, call ADLAPLACE_DEFINE_BACKEND and return
// parameters_ad_pack_ptr() or observations_ad_pack_ptr() from an export.

#include "adlaplace/extension.hpp"
#include "adlaplace/eval_impl.hpp"
#include "adlaplace/register_impl.hpp"

#endif
