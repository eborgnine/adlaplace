#ifndef HANDLE_CREATE_HPP
#define HANDLE_CREATE_HPP

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



static adlaplace_adpack_handle* get_adpack_handle_cpp(
	std::vector<GroupPack>& adPack,
		const HessianPack &hessian_inner,
	const HessianPack &hessian_outer,
  const Config& config)
{
  BackendContext* ctx = getBackendContext(&adPack, hessian_inner, hessian_outer, config);

  auto* h = new adlaplace_adpack_handle;
  h->api = &AD_API;
  h->ctx = static_cast<void*>(ctx);

  return h;
}


static void handle_finalizer(SEXP ext) {
	adlaplace_adpack_handle* h = (adlaplace_adpack_handle*)R_ExternalPtrAddr(ext);
	if (!h) return;
	if (h->api && h->api->destroy && h->ctx) h->api->destroy(h->ctx);
	delete h;
	R_ClearExternalPtr(ext);
}

static void backend_destroy(void* vctx) {
	BackendContext* ctx = (BackendContext*)vctx;
	delete ctx;
}






#endif