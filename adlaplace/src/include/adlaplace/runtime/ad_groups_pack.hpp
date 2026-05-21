#ifndef ADLAPLACE_AD_GROUPS_PACK_HPP
#define ADLAPLACE_AD_GROUPS_PACK_HPP

#include <Rcpp.h>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/creators/rviews.hpp"

inline hessian_template hessian_template_from_dgc(
  const DgCView& tpl,
  const char* label) {

  if (tpl.p.size() == 0) {
    Rcpp::stop("%s Hessian template missing; attach hessian_map result first", label);
  }

  const int n = tpl.nrow();
  const R_xlen_t p_len = tpl.p.size();
  const R_xlen_t i_len = tpl.i.size();
  if (static_cast<R_xlen_t>(n + 1) != p_len) {
    Rcpp::stop(
      "%s Hessian p_len=%d but expected nrow+1=%d",
      label,
      static_cast<int>(p_len),
      n + 1
    );
  }

  std::vector<int> p(static_cast<size_t>(p_len));
  std::vector<int> i(static_cast<size_t>(i_len));
  std::vector<int> x(static_cast<size_t>(i_len));
  for (R_xlen_t k = 0; k < p_len; ++k) {
    p[static_cast<size_t>(k)] = tpl.p[k];
  }
  for (R_xlen_t k = 0; k < i_len; ++k) {
    i[static_cast<size_t>(k)] = tpl.i[k];
    x[static_cast<size_t>(k)] = static_cast<int>(tpl.value<double>(k));
  }

  Eigen::Map<const hessian_template> upper(
    n, n, static_cast<int>(i_len),
    p.data(), i.data(), x.data()
  );
  hessian_template out = upper.selfadjointView<Eigen::Upper>();
  out.makeCompressed();
  return out;
}

// Copy Hessian sparsity templates and maps from a get_ad_fun()-style list
// (components outer, inner, map_outer, map_inner, sizes from hessian_map()).
inline void ad_groups_attach_hessians_from_list(
  ad_groups& groups,
  const Rcpp::List& ad_fun)
{
  if (groups.hessians_attached) {
    return;
  }

  if (!ad_fun.containsElementNamed("outer") || !ad_fun.containsElementNamed("inner")) {
    Rcpp::stop("ad_fun list must contain 'outer' and 'inner' (from hessian_map())");
  }
  if (!ad_fun.containsElementNamed("map_outer") || !ad_fun.containsElementNamed("map_inner")) {
    Rcpp::stop("ad_fun list must contain 'map_outer' and 'map_inner'");
  }
  if (!ad_fun.containsElementNamed("sizes")) {
    Rcpp::stop("ad_fun list must contain 'sizes'");
  }

  // Parameter block sizes (beta, gamma, theta); view into the list's sizes vector.
  groups.sizes = IntVecView(static_cast<SEXP>(ad_fun["sizes"]));
  if (!groups.sizes.has_name("beta") ||
      !groups.sizes.has_name("gamma") ||
      !groups.sizes.has_name("theta")) {
    Rcpp::stop("sizes must be a named vector: beta, gamma, theta");
  }

  // Outer (full) Hessian sparsity template -> Eigen int matrix on ad_groups.
  groups.hessian_outer = hessian_template_from_dgc(
    DgCView(Rcpp::as<Rcpp::S4>(ad_fun["outer"])), "outer");

  // Inner (gamma block) Hessian sparsity template -> Eigen int matrix on ad_groups.
  groups.hessian_inner = hessian_template_from_dgc(
    DgCView(Rcpp::as<Rcpp::S4>(ad_fun["inner"])), "inner");

  // Per-shard local-to-global index maps for outer and inner Hessian accumulation.
  groups.map_outer = hessian_map_view(Rcpp::as<Rcpp::List>(ad_fun["map_outer"]));
  groups.map_inner = hessian_map_view(Rcpp::as<Rcpp::List>(ad_fun["map_inner"]));

  groups.hessians_attached = true;
}

#endif
