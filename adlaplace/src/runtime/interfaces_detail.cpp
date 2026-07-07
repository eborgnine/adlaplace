#include "adlaplace/runtime/interfaces_detail.hpp"
#include "adlaplace/linalg/chol_update.hpp"

namespace {

static void validate_ad_fun_eval(const ad_fun& backend) {
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

adlaplace_adpack_handle* shard_handle(ad_fun* backend, size_t shard) {
  if (!backend) Rcpp::stop("ad_fun handle is NULL");
  if (shard >= backend->fun.size()) {
    Rcpp::stop("shard index %d out of range [0, %d]", (int)shard, (int)backend->fun.size() - 1);
  }
  adlaplace_adpack_handle* h = backend->fun[shard];
  if (!h) Rcpp::stop("ad_fun.fun[%d] is NULL", (int)shard);
  if (!h->api) Rcpp::stop("ad_fun.fun[%d] has NULL api", (int)shard);
  if (!h->ctx) Rcpp::stop("ad_fun.fun[%d] has NULL ctx", (int)shard);
  if (!h->api->f) Rcpp::stop("ad_fun.fun[%d] api->f is NULL", (int)shard);
  return h;
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

Rcpp::List sparsity_shard_from_handle(adlaplace_adpack_handle* h) {
  if (!h->api->get_sizes) {
    Rcpp::stop("backend api->get_sizes is NULL");
  }
  if (!h->api->get_sparse_pattern) {
    Rcpp::stop("backend api->get_sparse_pattern is NULL");
  }

  int n_inner = 0;
  int n_outer = 0;
  int n_beta = 0;
  int n_theta = 0;
  int nnz_grad_inner = 0;
  int nnz_grad_outer = 0;
  int nnz_hes_inner = 0;
  int nnz_hes_outer = 0;
  if (h->api->get_sizes(
        h->ctx, &n_inner, &n_outer, &n_beta, &n_theta,
        &nnz_grad_inner, &nnz_grad_outer,
        &nnz_hes_inner, &nnz_hes_outer) != 0) {
    Rcpp::stop("backend api->get_sizes failed");
  }

  std::vector<int> grad_inner(static_cast<size_t>(nnz_grad_inner));
  std::vector<int> grad_outer(static_cast<size_t>(nnz_grad_outer));
  std::vector<int> row_hes_inner(static_cast<size_t>(nnz_hes_inner));
  std::vector<int> col_hes_inner(static_cast<size_t>(nnz_hes_inner));
  std::vector<int> row_hes_outer(static_cast<size_t>(nnz_hes_outer));
  std::vector<int> col_hes_outer(static_cast<size_t>(nnz_hes_outer));

  if (h->api->get_sparse_pattern(
        h->ctx,
        grad_inner.data(), grad_outer.data(),
        row_hes_inner.data(), col_hes_inner.data(),
        row_hes_outer.data(), col_hes_outer.data()) != 0) {
    Rcpp::stop("backend api->get_sparse_pattern failed");
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
