#include "adlaplace/runtime.hpp"
#include "adlaplace/register.hpp"
#include "adlaplace/chol_update.hpp"

extern adlaplace_shard* adlaplace_make_shard(GroupPack&&);

ad_fun* clone_ad_fun(const ad_fun* src) {
  if (!src) {
    Rcpp::stop("clone_ad_fun: source handle is NULL");
  }
  if (src->fun.empty()) {
    Rcpp::stop("clone_ad_fun: source has no AD shards");
  }

  std::size_t n_beta = 0;
  std::size_t n_theta = 0;
  std::vector<GroupPack> packs;
  std::vector<ShardFactory> factories;
  packs.reserve(src->fun.size());
  factories.reserve(src->fun.size());

  for (adlaplace_shard* shard : src->fun) {
    if (!shard) {
      continue;
    }
    if (packs.empty()) {
      n_beta = shard->pack.n_beta;
      n_theta = shard->pack.n_theta;
    }
    factories.push_back(shard->factory ? shard->factory : adlaplace_make_shard);
    packs.push_back(clone_group_pack(shard->pack));
  }

  if (packs.empty()) {
    Rcpp::stop("clone_ad_fun: source has no valid shards");
  }

  ad_fun* copy = new ad_fun();
  copy->abi_version = src->abi_version;
  copy->fun.reserve(packs.size());
  for (size_t g = 0; g < packs.size(); ++g) {
    copy->fun.push_back(factories[g](std::move(packs[g])));
  }
  copy->configured_num_threads = src->configured_num_threads;
  copy->num_threads_configured = src->num_threads_configured;
  return copy;
}

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

namespace {

static void check_ad_fun_abi(const ad_fun& backend) {
  if (backend.abi_version != ADLAPLACE_ABI_VERSION) {
    Rcpp::stop(
      "ad_fun ABI mismatch (got %d, need %d): rebuild the extension package against this adlaplace",
      backend.abi_version,
      ADLAPLACE_ABI_VERSION
    );
  }
}

static void validate_ad_fun_eval(const ad_fun& backend) {
  check_ad_fun_abi(backend);
  if (backend.fun.empty()) {
    Rcpp::stop("ad_fun has no AD shards");
  }
}

static void validate_ad_fun_laplace(const ad_fun& backend) {
  validate_ad_fun_eval(backend);
  if (backend.sizes.sexp == R_NilValue ||
      !backend.sizes.has_name("beta") ||
      !backend.sizes.has_name("gamma") ||
      !backend.sizes.has_name("theta")) {
    Rcpp::stop("ad_fun missing sizes; call ad_fun() before inner_opt()");
  }
  if (!backend.hessians_attached) {
    Rcpp::stop("ad_fun missing Hessian templates; call ad_fun() before inner_opt()");
  }
}

static Rcpp::List ad_fun_list_from_s4(const Rcpp::S4& obj) {
  return Rcpp::List::create(
    Rcpp::Named("ad_fun") = obj.slot("ptr"),
    Rcpp::Named("group_sparsity") = obj.slot("group_sparsity"),
    Rcpp::Named("outer") = obj.slot("outer"),
    Rcpp::Named("inner") = obj.slot("inner"),
    Rcpp::Named("map_outer") = obj.slot("map_outer"),
    Rcpp::Named("map_inner") = obj.slot("map_inner"),
    Rcpp::Named("sizes") = obj.slot("sizes"),
    Rcpp::Named("chol_inner_list") = obj.slot("chol_inner_list")
  );
}

}  // namespace

adlaplace_shard* shard_handle(ad_fun* backend, size_t shard) {
  if (!backend) Rcpp::stop("ad_fun handle is NULL");
  if (shard >= backend->fun.size()) {
    Rcpp::stop("shard index %d out of range [0, %d]", (int)shard, (int)backend->fun.size() - 1);
  }
  adlaplace_shard* s = backend->fun[shard];
  if (!s) Rcpp::stop("ad_fun.fun[%d] is NULL", (int)shard);
  return s;
}

SEXP ad_fun_handle_sexp(const Rcpp::List& ad_fun_list) {
  if (ad_fun_list.containsElementNamed("ad_fun")) {
    return ad_fun_list["ad_fun"];
  }
  Rcpp::stop("ad_fun list must contain component 'ad_fun'");
}

ad_fun* ad_fun_from_list(const Rcpp::List& ad_fun_list) {
  SEXP handle_ptr = ad_fun_handle_sexp(ad_fun_list);
  auto* backend = static_cast<::ad_fun*>(R_ExternalPtrAddr(handle_ptr));
  if (!backend) {
    Rcpp::stop("ad_fun external pointer is NULL (cleared?)");
  }

  ad_fun_attach_hessians_from_list(*backend, ad_fun_list);
  return backend;
}

ad_fun* ad_fun_from_handle(SEXP handle) {
  auto* backend = static_cast<::ad_fun*>(R_ExternalPtrAddr(handle));
  if (!backend) {
    Rcpp::stop("ad_fun external pointer is NULL (cleared?)");
  }
  check_ad_fun_abi(*backend);
  return backend;
}

ad_fun* resolve_ad_fun_eval(SEXP ad_fun_ptr) {
  if (Rf_isNull(ad_fun_ptr)) {
    Rcpp::stop("ad_fun_ptr must not be NULL");
  }
  if (!Rf_inherits(ad_fun_ptr, "ad_fun_ptr")) {
    Rcpp::stop("expected ad_fun_ptr");
  }
  ad_fun* backend = ad_fun_from_handle(ad_fun_ptr);
  validate_ad_fun_eval(*backend);
  return backend;
}

ad_fun* resolve_ad_fun_laplace(const Rcpp::S4& ad_fun_s4) {
  if (!ad_fun_s4.inherits("ad_fun")) {
    Rcpp::stop("expected an ad_fun object; call ad_fun(ad_fun_ptr) first");
  }
  ::ad_fun* backend = ad_fun_from_list(ad_fun_list_from_s4(ad_fun_s4));
  validate_ad_fun_laplace(*backend);
  return backend;
}

std::vector<size_t> resolve_shard_indices(
  size_t n_shards,
  const Rcpp::IntegerVector& shards) {
  std::vector<size_t> shard_idx;
  if (shards.size() == 0) {
    shard_idx.resize(n_shards);
    for (size_t s = 0; s < n_shards; ++s) shard_idx[s] = s;
    return shard_idx;
  }

  shard_idx.reserve(shards.size());
  for (R_xlen_t k = 0; k < shards.size(); ++k) {
    if (shards[k] == NA_INTEGER) {
      Rcpp::stop("shards contains NA at position %d", (int)k + 1);
    }
    if (shards[k] < 0 || static_cast<size_t>(shards[k]) >= n_shards) {
      Rcpp::stop("shards index %d out of range [0, %d]", shards[k], (int)n_shards - 1);
    }
    shard_idx.push_back(static_cast<size_t>(shards[k]));
  }
  return shard_idx;
}

Rcpp::List sparsity_shard_from_handle(adlaplace_shard* shard) {
  int n_inner = 0;
  int n_outer = 0;
  int n_beta = 0;
  int n_theta = 0;
  int nnz_grad_inner = 0;
  int nnz_grad_outer = 0;
  int nnz_hes_inner = 0;
  int nnz_hes_outer = 0;
  if (shard->get_sizes(
        &n_inner, &n_outer, &n_beta, &n_theta,
        &nnz_grad_inner, &nnz_grad_outer,
        &nnz_hes_inner, &nnz_hes_outer) != 0) {
    Rcpp::stop("shard get_sizes failed");
  }

  std::vector<int> grad_inner(static_cast<size_t>(nnz_grad_inner));
  std::vector<int> grad_outer(static_cast<size_t>(nnz_grad_outer));
  std::vector<int> row_hes_inner(static_cast<size_t>(nnz_hes_inner));
  std::vector<int> col_hes_inner(static_cast<size_t>(nnz_hes_inner));
  std::vector<int> row_hes_outer(static_cast<size_t>(nnz_hes_outer));
  std::vector<int> col_hes_outer(static_cast<size_t>(nnz_hes_outer));

  if (shard->get_sparse_pattern(
        grad_inner.data(), grad_outer.data(),
        row_hes_inner.data(), col_hes_inner.data(),
        row_hes_outer.data(), col_hes_outer.data()) != 0) {
    Rcpp::stop("shard get_sparse_pattern failed");
  }

  return Rcpp::List::create(
    Rcpp::Named("grad") = Rcpp::IntegerVector(
      grad_outer.begin(), grad_outer.end()),
    Rcpp::Named("grad_inner") = Rcpp::IntegerVector(
      grad_inner.begin(), grad_inner.end()),
    Rcpp::Named("row_hess") = Rcpp::IntegerVector(
      row_hes_outer.begin(), row_hes_outer.end()),
    Rcpp::Named("col_hess") = Rcpp::IntegerVector(
      col_hes_outer.begin(), col_hes_outer.end()),
    Rcpp::Named("row_hess_inner") = Rcpp::IntegerVector(
      row_hes_inner.begin(), row_hes_inner.end()),
    Rcpp::Named("col_hess_inner") = Rcpp::IntegerVector(
      col_hes_inner.begin(), col_hes_inner.end())
  );
}
