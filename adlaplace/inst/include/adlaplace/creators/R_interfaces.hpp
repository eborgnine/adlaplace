#ifdef R_INTERFACES_HPP
#define R_INTERFACES_HPP

#include<Rcpp.h>

#include "adlaplace/adlaplace.hpp" 

static const adlaplace_adpack_api AD_API = {
	ADLAPLACE_ADPACK_API_VERSION,
  1,                // thread_safe
  0,                // Ngamma (set at runtime in handle if you prefer)
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
	if (h->api && h->api->destroy && h->ctx) h->api->destroy(h->ctx);
	delete h;
	R_ClearExternalPtr(ext);
}

// helper: get handle safely
static inline adlaplace_adpack_handle* get_handle(SEXP handle_sexp) {
	auto* h = static_cast<adlaplace_adpack_handle*>(R_ExternalPtrAddr(handle_sexp));
	if (!h) Rcpp::stop("backendContext handle is NULL (external pointer cleared?)");
	if (!h->api) Rcpp::stop("backendContext handle has NULL api");
	if (!h->ctx) Rcpp::stop("backendContext handle has NULL ctx");
	if (!h->api->f) Rcpp::stop("backendContext api->f is NULL");
	return h;
}

SEXP getAdFun_h(
	const Rcpp::List data, 
	const Rcpp::List config) {

	const Data dataC(data);
	const Config configC(config);

	std::vector<GroupPack> adPack = getAdFun(dataC, configC);
	const Rcpp::List sparsity_list = extractSparsity(adPack);

	Rcpp::Environment ns = Rcpp::Environment::namespace_env("adlaplace");
	Rcpp::Function hessianMap = ns["hessianMap"];
	const Rcpp::List hessians = Rcpp::as<Rcpp::List>(hessianMap(sparsity_list, config));
	
	std::array<HessianPack,2> hessiansC = hessianPackFromList(hessians);

	BackendContext* ctx = getBackendContext(&adPack, hessianC[0], hessianC[1], configC);

	auto* h = new adlaplace_adpack_handle;
	h->api = &AD_API;
	h->ctx = static_cast<void*>(ctx);

	SEXP handle2 = R_MakeExternalPtr((void*)h, R_NilValue);
	R_RegisterCFinalizerEx(handle2, handle_finalizer, TRUE);
	Rf_setAttrib(handle2, R_ClassSymbol, Rf_mkString("adlaplace_handle_ptr"));

	return handle2;
}

double jointLogDens_h(Rcpp::NumericVector x, SEXP handle) {

	adlaplace_adpack_handle* h = get_handle(backendContext);

	size_t Nparams, Ngroups, Nbeta, Ngamma, Ntheta;  
	h->get_sizes(&Nparams, &Ngroups, &Nbeta, & Ngamma, & Ntheta);

	const size_t Ngroups = h->Ngroups;
	const size_t Nparams = h->Nparams;

	if (static_cast<size_t>(x.size()) != Nparams) {
		Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
	}

  // Evaluate each group's contribution
	for (size_t g = 0; g < Ngroups; ++g) {
		double fg = 0.0;
		int gi = static_cast<int>(g);
		int rc = h->api->f(h->ctx, &gi, x.begin(), &fg);
		if (rc != 0) {
			Rcpp::stop("backend api->f failed for group %d with code %d", gi, rc);
		}
	}
	return fg;
}
#endif