#ifndef ADLAPLACE_AD_GROUPS_PACK_HPP
#define ADLAPLACE_AD_GROUPS_PACK_HPP

#include <Rcpp.h>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/creators/rviews.hpp"

hessian_template hessian_template_from_dgc(
  const DgCView& tpl,
  const char* label);

void ad_fun_attach_hessians_from_list(
  ad_fun& shards,
  const Rcpp::List& ad_fun);

#endif
