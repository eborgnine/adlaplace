#ifndef ADLAPLACE_R_INTERFACES_HPP
#define ADLAPLACE_R_INTERFACES_HPP

#include <Rcpp.h>
#include <memory>

#include "adlaplace/adlaplace.hpp"
#include "adlaplace/eval/fgh.hpp"
#include "adlaplace/creators/handle.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"
#include "adlaplace/api/stubs.hpp"

std::vector<GroupPack> getAdFun(const Data& data, const Config& config);
Rcpp::List extract_sparsity(const std::vector<GroupPack> &adFun);

static const adlaplace_adpack_api AD_API = {
	ADLAPLACE_ADPACK_API_VERSION,
	1,
	&eval_f,
	&eval_grad,
	&eval_hess,
	&get_sizes,
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


inline Rcpp::List getAdFun_h(
	const Rcpp::List &data,
	const Rcpp::List &config) {

	const Data dataC(data);
	const Config configC(config);

	auto adFun = std::make_unique<std::vector<GroupPack>>(getAdFun(dataC, configC));
	const Rcpp::List sparsity_list = extract_sparsity(*adFun);

	Rcpp::Environment ns = Rcpp::Environment::namespace_env("adlaplace");
	Rcpp::Function hessianMap = ns["hessianMap"];
	const Rcpp::List hessians = Rcpp::as<Rcpp::List>(hessianMap(sparsity_list, config));

	std::array<HessianPack, 2> hessiansC = adlaplace_hessianPackFromList(hessians);
	auto ctx = std::unique_ptr<BackendContext>(
		get_backend_context(adFun.get(), hessiansC[0], hessiansC[1], configC)
	);

	auto h = std::unique_ptr<adlaplace_adpack_handle>(new adlaplace_adpack_handle);
	h->api = &AD_API;
	h->ctx = static_cast<void*>(ctx.get());

	SEXP handle = R_MakeExternalPtr((void*)h.get(), R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(handle, handle_finalizer, TRUE);
	Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));

	// ownership moves to finalizer path
	adFun.release();
	ctx.release();
	h.release();
	Rcpp::List result = Rcpp::List::create(
		Rcpp::Named("adFun") = handle,
		Rcpp::Named("sparsity") = sparsity_list,
		Rcpp::Named("hessians") = hessians
		);

	return result;
}

#endif
