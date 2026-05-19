#ifndef ADLAPLACE_R_INTERFACES_HPP
#define ADLAPLACE_R_INTERFACES_HPP

#include <Rcpp.h>
#include <memory>

#include "adlaplace/adlaplace.hpp"
#include "adlaplace/eval/fgh.hpp"
#include "adlaplace/creators/handle.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"

std::vector<GroupPack> getAdFun(const Data& data, const Config& config);
Rcpp::List extract_sparsity(const std::vector<GroupPack> &adFun);

static const adlaplace_adpack_api AD_API = {
	ADLAPLACE_ADPACK_API_VERSION,
	1,
	&eval_f,
	&eval_grad,
	&eval_hess,
	&get_sizes,
	&get_sparse_sizes,
	&get_sparse_pattern,
	&get_hessian,
	&eval_trace_hinv_t,
	&backend_destroy,
	NULL
};

static void handle_finalizer(SEXP ext) {
	adlaplace_adpack_handle* h = (adlaplace_adpack_handle*)R_ExternalPtrAddr(ext);
	if (!h) return;
	if (h->api && h->api->destroy && h->ctx) {
		h->api->destroy(h->ctx);
	}
	delete h;
	R_ClearExternalPtr(ext);
}

static inline adlaplace_adpack_handle* get_handle(SEXP handle_sexp) {
	SEXP handle_ptr = handle_sexp;
	if (TYPEOF(handle_sexp) == VECSXP) {
		Rcpp::List maybe_list(handle_sexp);
		if (maybe_list.containsElementNamed("adFun")) {
			handle_ptr = maybe_list["adFun"];
		}
	}

	auto* h = static_cast<adlaplace_adpack_handle*>(R_ExternalPtrAddr(handle_ptr));
	if (!h) Rcpp::stop("backendContext handle is NULL (external pointer cleared?)");
	if (!h->api) Rcpp::stop("backendContext handle has NULL api");
	if (!h->ctx) Rcpp::stop("backendContext handle has NULL ctx");
	if (!h->api->f) Rcpp::stop("backendContext api->f is NULL");
	return h;
}

static inline std::vector<size_t> resolve_groups(
	size_t Ngroups,
	const Rcpp::IntegerVector& Sgroups) {
	std::vector<size_t> groups;
	if (Sgroups.size() == 0) {
		groups.resize(Ngroups);
		for (size_t g = 0; g < Ngroups; ++g) groups[g] = g;
		return groups;
	}

	groups.reserve(Sgroups.size());
	for (R_xlen_t k = 0; k < Sgroups.size(); ++k) {
		if (Sgroups[k] == NA_INTEGER) {
			Rcpp::stop("Sgroups contains NA at position %d", (int)k + 1);
		}
		if (Sgroups[k] < 0 || static_cast<size_t>(Sgroups[k]) >= Ngroups) {
			Rcpp::stop("Sgroups index %d out of range [0, %d]", Sgroups[k], (int)Ngroups - 1);
		}
		groups.push_back(static_cast<size_t>(Sgroups[k]));
	}
	return groups;
}

inline Rcpp::List sparsity_shard_from_handle(adlaplace_adpack_handle* h, int group) {
	if (!h->api->get_sparse_sizes) {
		Rcpp::stop("backend api->get_sparse_sizes is NULL");
	}
	if (!h->api->get_sparse_pattern) {
		Rcpp::stop("backend api->get_sparse_pattern is NULL");
	}

	int gi = group;
	int n_inner = 0;
	int n_outer = 0;
	int nnz_grad_inner = 0;
	int nnz_grad_outer = 0;
	int nnz_hes_inner = 0;
	int nnz_hes_outer = 0;
	if (h->api->get_sparse_sizes(
			h->ctx, &gi, &n_inner, &n_outer,
			&nnz_grad_inner, &nnz_grad_outer,
			&nnz_hes_inner, &nnz_hes_outer) != 0) {
		Rcpp::stop("backend api->get_sparse_sizes failed for group %d", gi);
	}

	std::vector<int> grad_inner(static_cast<size_t>(nnz_grad_inner));
	std::vector<int> grad_outer(static_cast<size_t>(nnz_grad_outer));
	std::vector<int> row_hes_inner(static_cast<size_t>(nnz_hes_inner));
	std::vector<int> col_hes_inner(static_cast<size_t>(nnz_hes_inner));
	std::vector<int> row_hes_outer(static_cast<size_t>(nnz_hes_outer));
	std::vector<int> col_hes_outer(static_cast<size_t>(nnz_hes_outer));

	if (h->api->get_sparse_pattern(
			h->ctx, &gi,
			grad_inner.data(), grad_outer.data(),
			row_hes_inner.data(), col_hes_inner.data(),
			row_hes_outer.data(), col_hes_outer.data()) != 0) {
		Rcpp::stop("backend api->get_sparse_pattern failed for group %d", gi);
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

inline Rcpp::List sparsity_list_from_handle(adlaplace_adpack_handle* h) {
	AdGroups* groups = groups_ctx(h->ctx);
	const size_t Ngroups = groups->size();
	Rcpp::List sparsity(Ngroups);
	for (size_t g = 0; g < Ngroups; ++g) {
		sparsity[static_cast<R_xlen_t>(g)] = sparsity_shard_from_handle(h, static_cast<int>(g));
	}
	return sparsity;
}

inline SEXP buildAdGroups_h(
	const Rcpp::List& data,
	const Rcpp::List& config) {

	const Data dataC(data);
	const Config configC(config);

	auto* groups = new AdGroups(getAdFun(dataC, configC));

	auto* h = new adlaplace_adpack_handle();
	h->api = &AD_API;
	h->ctx = static_cast<void*>(groups);

	SEXP handle = R_MakeExternalPtr(static_cast<void*>(h), R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(handle, handle_finalizer, TRUE);
	Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));
	return handle;
}

inline void finalizeAdHandle_h(
	SEXP handle_sexp,
	const Rcpp::List& hessians) {
	(void)handle_sexp;
	(void)hessians;
}

inline Rcpp::List get_sizes_from_handle(adlaplace_adpack_handle* h) {
	size_t Nparams = 0;
	size_t Ngroups = 0;
	size_t Nbeta = 0;
	size_t Ngamma = 0;
	size_t Ntheta = 0;
	if (!h->api->get_sizes) {
		Rcpp::stop("backend api->get_sizes is NULL");
	}
	if (h->api->get_sizes(h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta) != 0) {
		Rcpp::stop("backend api->get_sizes failed");
	}
	return Rcpp::List::create(
		Rcpp::Named("Nparams") = static_cast<int>(Nparams),
		Rcpp::Named("Ngroups") = static_cast<int>(Ngroups)
	);
}

inline SEXP getAdFun_h(
	const Rcpp::List &data,
	const Rcpp::List &config) {
	return buildAdGroups_h(data, config);
}

#endif
