#include "adlaplace/runtime/ad_fun_pack.hpp"
#include "adlaplace/linalg/chol_update.hpp"

hessian_template hessian_template_from_dgc(
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

void ad_fun_attach_hessians_from_list(
  ad_fun& shards,
  const Rcpp::List& ad_fun) {

  if (!ad_fun.containsElementNamed("outer") || !ad_fun.containsElementNamed("inner")) {
    Rcpp::stop("ad_fun list must contain 'outer' and 'inner' (from hessian_map())");
  }
  if (!ad_fun.containsElementNamed("map_outer") || !ad_fun.containsElementNamed("map_inner")) {
    Rcpp::stop("ad_fun list must contain 'map_outer' and 'map_inner'");
  }
  if (!ad_fun.containsElementNamed("sizes")) {
    Rcpp::stop("ad_fun list must contain 'sizes'");
  }

  shards.sizes = IntVecView(static_cast<SEXP>(ad_fun["sizes"]));
  if (!shards.sizes.has_name("beta") ||
      !shards.sizes.has_name("gamma") ||
      !shards.sizes.has_name("theta")) {
    Rcpp::stop("sizes must be a named vector: beta, gamma, theta");
  }

  if (shards.hessians_attached) {
    return;
  }

  shards.hessian_outer = hessian_template_from_dgc(
    DgCView(Rcpp::as<Rcpp::S4>(ad_fun["outer"])), "outer");

  shards.hessian_inner = hessian_template_from_dgc(
    DgCView(Rcpp::as<Rcpp::S4>(ad_fun["inner"])), "inner");

  shards.map_outer = hessian_map_view(Rcpp::as<Rcpp::List>(ad_fun["map_outer"]));
  shards.map_inner = hessian_map_view(Rcpp::as<Rcpp::List>(ad_fun["map_inner"]));

  shards.hessians_attached = true;
  ad_fun_attach_chol_pattern_from_list(shards, ad_fun);
}
