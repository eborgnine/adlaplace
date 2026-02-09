#ifndef ADLAPLACE_R_INTERFACES_HPP
#define ADLAPLACE_R_INTERFACES_HPP

#include <Rcpp.h>
#include <memory>

#include "adlaplace/adlaplace.hpp"
#include "adlaplace/creators/api_functions.hpp"
#include "adlaplace/creators/handle.hpp"

std::vector<GroupPack> getAdFun(const Data& data, const Config& config);
Rcpp::List extract_sparsity(const std::vector<GroupPack> &adPack);

static const adlaplace_adpack_api AD_API = {
	ADLAPLACE_ADPACK_API_VERSION,
	1,
	&eval_f,
	&eval_grad,
	&eval_hess,
	&get_sizes,
	&get_hessian,
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

	auto adPack = std::make_unique<std::vector<GroupPack>>(getAdFun(dataC, configC));
	const Rcpp::List sparsity_list = extract_sparsity(*adPack);

	Rcpp::Environment ns = Rcpp::Environment::namespace_env("adlaplace");
	Rcpp::Function hessianMap = ns["hessianMap"];
	const Rcpp::List hessians = Rcpp::as<Rcpp::List>(hessianMap(sparsity_list, config));

	std::array<HessianPack, 2> hessiansC = hessianPackFromList(hessians);
	auto ctx = std::unique_ptr<BackendContext>(
		get_backend_context(adPack.get(), hessiansC[0], hessiansC[1], configC)
	);

	auto h = std::unique_ptr<adlaplace_adpack_handle>(new adlaplace_adpack_handle);
	h->api = &AD_API;
	h->ctx = static_cast<void*>(ctx.get());

	SEXP handle = R_MakeExternalPtr((void*)h.get(), R_NilValue, R_NilValue);
	R_RegisterCFinalizerEx(handle, handle_finalizer, TRUE);
	Rf_setAttrib(handle, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));

	// ownership moves to finalizer path
	adPack.release();
	ctx.release();
	h.release();
	Rcpp::List result = Rcpp::List::create(
		Rcpp::Named("adFun") = handle,
		Rcpp::Named("sparsity") = sparsity_list,
		Rcpp::Named("hessians") = hessians
		);

	return result;
}

inline double jointLogDens_h(
	const Rcpp::NumericVector &x,
	SEXP handle,
	SEXP Sgroups = R_NilValue) {

	adlaplace_adpack_handle* h = get_handle(handle);

	size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
	const int rc_sizes = h->api->get_sizes(
		h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
	);
	if (rc_sizes != 0) {
		Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
	}

	if (static_cast<size_t>(x.size()) != Nparams) {
		Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
	}

	const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
		? Rcpp::IntegerVector()
		: Rcpp::as<Rcpp::IntegerVector>(Sgroups);
	const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);
	double total = 0.0;
	for (size_t g : groups) {

		double fg = 0.0;
		int gi = static_cast<int>(g);
		const int rc = h->api->f(h->ctx, &gi, x.begin(), &fg);
		if (rc != 0) {
			Rcpp::stop("backend api->f failed for group %d with code %d", gi, rc);
		}
		total += fg;
	}

	return total;
}

inline Rcpp::NumericVector grad_h(
	const Rcpp::NumericVector &x,
	SEXP handle,
	const bool inner = false,
	SEXP Sgroups = R_NilValue) {

	adlaplace_adpack_handle* h = get_handle(handle);
	if (!h->api->f_grad) {
		Rcpp::stop("backendContext api->f_grad is NULL");
	}

	size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
	const int rc_sizes = h->api->get_sizes(
		h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
	);
	if (rc_sizes != 0) {
		Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
	}

	if (static_cast<size_t>(x.size()) != Nparams) {
		Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
	}

	Rcpp::NumericVector grad_out(Nparams, 0.0);
	double f_total = 0.0;
	const bool inner_local = inner;
	const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
		? Rcpp::IntegerVector()
		: Rcpp::as<Rcpp::IntegerVector>(Sgroups);
	const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

	for (size_t g : groups) {
		int gi = static_cast<int>(g);
		const int rc = h->api->f_grad(
			h->ctx, &gi, x.begin(), &inner_local, &f_total, grad_out.begin()
		);
		if (rc != 0) {
			Rcpp::stop("backend api->f_grad failed for group %d with code %d", gi, rc);
		}
	}

	return grad_out;
}

inline Rcpp::S4 hess_h(
	const Rcpp::NumericVector &x,
	SEXP handle,
	const bool inner = false,
	SEXP Sgroups = R_NilValue) {

	adlaplace_adpack_handle* h = get_handle(handle);
	if (!h->api->f_grad_hess) {
		Rcpp::stop("backendContext api->f_grad_hess is NULL");
	}
	if (!h->api->get_hessian) {
		Rcpp::stop("backendContext api->get_hessian is NULL");
	}

	size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
	const int rc_sizes = h->api->get_sizes(
		h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
	);
	if (rc_sizes != 0) {
		Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
	}

	if (static_cast<size_t>(x.size()) != Nparams) {
		Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
	}

	const bool inner_local = inner;
	const int* hess_p = NULL;
	const int* hess_i = NULL;
	size_t hess_p_len = 0;
	size_t hess_i_len = 0;
	const int rc_hessian = h->api->get_hessian(
		h->ctx, &inner_local, &hess_p, &hess_p_len, &hess_i, &hess_i_len
	);
	if (rc_hessian != 0) {
		Rcpp::stop("backend api->get_hessian failed with code %d", rc_hessian);
	}

	Rcpp::NumericVector hess_out(hess_i_len, 0.0);
	Rcpp::NumericVector grad_scratch(Nparams, 0.0);
	double f_total = 0.0;
	const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
		? Rcpp::IntegerVector()
		: Rcpp::as<Rcpp::IntegerVector>(Sgroups);
	const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);


	for (size_t g : groups) {
		int gi = static_cast<int>(g);
		const int rc = h->api->f_grad_hess(
			h->ctx,
			&gi,
			x.begin(),
			&inner_local,
			&f_total,
			grad_scratch.begin(),
			hess_out.begin()
		);
		if (rc != 0) {
			Rcpp::stop("backend api->f_grad_hess failed for group %d with code %d", gi, rc);
		}
	}

	const int ncol = static_cast<int>(hess_p_len > 0 ? hess_p_len - 1 : 0);
	Rcpp::IntegerVector p_out(hess_p_len);
	Rcpp::IntegerVector i_out(hess_i_len);
	for (size_t k = 0; k < hess_p_len; ++k) p_out[k] = static_cast<int>(hess_p[k]);
	for (size_t k = 0; k < hess_i_len; ++k) i_out[k] = static_cast<int>(hess_i[k]);

	// Ensure Matrix namespace/classes are available before constructing S4 class.
	Rcpp::Environment matrix_ns = Rcpp::Environment::namespace_env("Matrix");
	(void)matrix_ns;
	Rcpp::S4 out("dsCMatrix");
	out.slot("i") = i_out;
	out.slot("p") = p_out;
	out.slot("x") = hess_out;
	out.slot("Dim") = Rcpp::IntegerVector::create(ncol, ncol);
	out.slot("uplo") = Rcpp::String("U");

	return out;
}

#endif
