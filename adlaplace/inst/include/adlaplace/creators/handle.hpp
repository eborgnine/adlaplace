#ifndef HANDLE_CREATE_HPP
#define HANDLE_CREATE_HPP

#include "adlaplace/runtime/backend.hpp"


static BackendContext* get_backend_context(
	std::vector<GroupPack> *adFun,
	const HessianPack &hessian_inner,
	const HessianPack &hessian_outer,
	const Config &config
	) {

  auto* ctx = new BackendContext;

	ctx->adFun = adFun;

	ctx->Nparams = config.Nparams;
	ctx->Ngroups = adFun ? adFun->size() : 0;
	ctx->Nbeta = config.Nbeta;
	ctx->Ngamma = config.Ngamma;
	ctx->Ntheta = config.Ntheta;

	ctx->hessian_inner = hessian_inner;

	ctx->hessian_outer = hessian_outer;

	return ctx;
}




static void backend_destroy(void* vctx) {
	BackendContext* ctx = (BackendContext*)vctx;
	if (!ctx) return;
	delete ctx->adFun;
	delete ctx;
}


#endif
