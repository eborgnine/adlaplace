#ifndef HANDLE_CREATE_HPP
#define HANDLE_CREATE_HPP

#include<adlaplace/backend.hpp>


static BackendContext* get_backend_context(
	std::vector<GroupPack> *adPack,
	const HessianPack &hessian_inner,
	const HessianPack &hessian_outer,
	const Config &config
	) {

  auto* ctx = new BackendContext;

	ctx->adPack = adPack;

	ctx->Nparams = config.Nparams;
	ctx->Ngroups = config.Ngroups;
	ctx->Nbeta = config.Nbeta;
	ctx->Ngamma = config.Ngamma;
	ctx->Ntheta = config.Ntheta;

	ctx->hessian_inner = hessian_inner;

	ctx->hessian_outer = hessian_outer;

	return ctx;
}




static void backend_destroy(void* vctx) {
	BackendContext* ctx = (BackendContext*)vctx;
	delete ctx;
}


#endif