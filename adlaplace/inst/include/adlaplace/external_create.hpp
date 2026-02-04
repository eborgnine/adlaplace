#ifndef EXTERNAL_CREATE_HPP
#define EXTERNAL_CREATE_HPP

#include <cppad/cppad.hpp>
#include<Rinternals.h>

#include<vector>
#include<cstddef>

#include "adlaplace/adpack_api.h"
#include "adlaplace/adpack_handle.h"


// backend private context
struct BackendContext {
	std::vector<GroupPack> *adPack;
	std::vector<std::vector<int>> map_inner; // p, local index, global index
	std::vector<std::vector<int>> map_outer; // 
};

static int eval_f(void* vctx, const int *i, const double* x, double* out_f) {
	auto* ctx = static_cast<BackendContext*>(vctx);
	if (*i < 0) return 2;
	size_t ist = (size_t)*i;

#ifdef DEBUG
	GroupPack& gp = ctx->adPack->at(ist);
#else	
	GroupPack &gp = (*(ctx->adPack))[ist];
#endif

	const size_t Nparams = gp.x.size();
	for(size_t D=0;D<Nparams;D++) {
		gp.x[D] = x[D];
	}	

	*out_f += gp.fun.Forward(0, gp.x)[0];
	return 0;
}

static int eval_grad_inner(void* vctx, const int *i, 
	const double* x, double* out_f, double* out_grad) {

	auto* ctx = static_cast<BackendContext*>(vctx);
	if (*i < 0) return 2;
	size_t ist = (size_t)*i;

#ifdef DEBUG
	GroupPack& gp = ctx->adPack->at(ist);
#else	
	GroupPack &gp = (*(ctx->adPack))[ist];
#endif


	const size_t Nparams = gp.x.size();
	for(size_t D=0;D<Nparams;D++) {
		gp.x[D] = x[D];
	}

	auto &pattern_here = gp.pattern_grad_inner;

	*out_f += gp.fun.Forward(0, gp.x)[0];
	gp.fun.sparse_jac_rev(
		gp.x,
		pattern_here,
		gp.unused_pattern,
		JAC_COLOR,
		gp.work_grad_inner);

	const size_t NoutGrad = pattern_here.nnz();
	const auto& cols = pattern_here.col();
	const auto& vals = pattern_here.val();
	for(size_t D=0; D < NoutGrad; D++) {
		out_grad[cols[D]] += vals[D];
	}

	return 0;
}

static int eval_hess_inner(void* vctx, const int *i, const double* x, double* out_f,
	double* out_grad, double* out_hess) {

	auto* ctx = static_cast<BackendContext*>(vctx);
	if (*i < 0) return 2;
	size_t ist = (size_t)*i;

#ifdef DEBUG
	GroupPack& gp = ctx->adPack->at(ist);
#else	
	GroupPack &gp = (*(ctx->adPack))[ist];
#endif

	const size_t Nparams = gp.x.size();
	for(size_t D=0;D<Nparams;D++) {
		gp.x[D] = x[D];
	}

	*out_f += gp.fun.Forward(0, gp.x)[0];

	auto &pattern_here_grad = gp.pattern_grad_inner;
	auto &pattern_here_hess = gp.pattern_hess_inner;
	const auto &csc_here = ctx->csc_inner;

	gp.fun.sparse_jac_rev(
		gp.x,
		pattern_here_grad,
		gp.unused_pattern,
		JAC_COLOR,
		gp.work_grad_inner);

	gp.fun.sparse_hes(
		gp.x,  
		gp.w,
		pattern_here_hess,              
		gp.unused_pattern,
		HESS_COLOR,                
		gp.work_hess_inner              
		);

	const size_t NoutGrad = pattern_here_grad.nnz();
	const auto& cols = pattern_here_grad.col();
	const auto& vals = pattern_here_grad.val();

	for(size_t D=0; D < NoutGrad; D++) {
		out_grad[cols[D]] += vals[D];
	}

	const std::size_t start = map_inner[0][ist];
	const std::size_t end   = map_inner[0][ist+1];;
	const std::size_t Nhere = end - start;
	const auto& vals_hess = pattern_here_hess.val();


	for(size_t Di=start;Di < end; Di++) {
		out_hess[map_inner[2][Di]] += vals_hess[map_inner[1][Di] ];
	}
	return 0;
}

static void backend_destroy(void* vctx) {
	BackendContext* ctx = (BackendContext*)vctx;
	delete ctx;
}

static const adlaplace_adpack_api AD_API = {
	ADLAPLACE_ADPACK_API_VERSION,
  1,                // thread_safe
  0,                // Ngamma (set at runtime in handle if you prefer)
  &eval_f,
  &eval_grad_inner,
  &eval_hess_inner,
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


#endif
