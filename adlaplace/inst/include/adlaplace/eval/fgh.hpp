#ifndef API_FUNCTIONS_HPP
#define API_FUNCTIONS_HPP

#include <cppad/cppad.hpp>
#include <Rinternals.h>

#include <vector>
#include <cstddef>

#include "adlaplace/runtime/backend.hpp"


static int get_sizes(void* vctx, size_t* Nparams, size_t* Ngroups,
	size_t* Nbeta, size_t* Ngamma, size_t* Ntheta){

	auto* ctx = static_cast<BackendContext*>(vctx);

	*Nparams = ctx->Nparams;
	*Ngroups = ctx->Ngroups;
	*Nbeta =   ctx->Nbeta;
	*Ngamma =  ctx->Ngamma;
	*Ntheta =  ctx->Ntheta;

	return 0;
}

static int get_hessian(void* vctx,
	const bool *inner,
	const int** p, size_t* p_len,
	const int** i, size_t* i_len){

	auto* ctx = static_cast<BackendContext*>(vctx);
	const bool innerv = *inner;

	*p = innerv ? ctx->hessian_inner.hessian_p.data() : ctx->hessian_outer.hessian_p.data();
	*i = innerv ? ctx->hessian_inner.hessian_i.data() : ctx->hessian_outer.hessian_i.data();
	*p_len = innerv ? ctx->hessian_inner.hessian_p.size() : ctx->hessian_outer.hessian_p.size();
	*i_len = innerv ? ctx->hessian_inner.hessian_i.size() : ctx->hessian_outer.hessian_i.size();

	return 0L;
}

static int eval_f(void* vctx, const int *i, const double* x, double* out_f) {
	auto* ctx = static_cast<BackendContext*>(vctx);
	if (*i < 0) return 2;
	size_t ist = (size_t)*i;
	if (ist >= ctx->adPack->size()) return 3;

	GroupPack &gp = (*(ctx->adPack))[ist];
	const size_t Nparams = gp.x.size();

	for(size_t D=0;D<Nparams;D++) {
		gp.x[D] = x[D];
	}	

	*out_f += gp.fun.Forward(0, gp.x)[0];
	return 0;
}

static int eval_grad(void* vctx, const int *i, 
	const double* x, const bool *inner, 
	double* out_f, double* out_grad) {

	auto* ctx = static_cast<BackendContext*>(vctx);
	if (*i < 0) return 2;
	size_t ist = (size_t)*i;
	if (ist >= ctx->adPack->size()) return 3;

	const bool innerv = *inner;

	GroupPack &gp = (*(ctx->adPack))[ist];


	const size_t Nparams = gp.x.size();
	for(size_t D=0;D<Nparams;D++) {
		gp.x[D] = x[D];
	}

	auto &pattern_here = innerv?gp.pattern_grad_inner:gp.pattern_grad;
	auto &work_here = innerv?gp.work_inner_grad:gp.work_grad;

	*out_f += gp.fun.Forward(0, gp.x)[0];
	gp.fun.sparse_jac_rev(
		gp.x,
		pattern_here,
		gp.unused_pattern,
		"cppad",
		work_here);

	const size_t NoutGrad = pattern_here.nnz();
	const auto& cols = pattern_here.col();
	const auto& vals = pattern_here.val();
	for(size_t D=0; D < NoutGrad; D++) {
		out_grad[cols[D]] += vals[D];
	}

	return 0;
}

static int eval_hess(void* vctx, const int *i, const double* x, 
	const bool *inner, double* out_f,
	double* out_grad, double* out_hess) {

	auto* ctx = static_cast<BackendContext*>(vctx);
	if (*i < 0) return 2;
	size_t ist = (size_t)*i;
	if (ist >= ctx->adPack->size()) return 3;

	const bool innerv = *inner;

	GroupPack &gp = (*(ctx->adPack))[ist];

	const size_t Nparams = gp.x.size();
	for(size_t D=0;D<Nparams;D++) {
		gp.x[D] = x[D];
	}

	*out_f += gp.fun.Forward(0, gp.x)[0];

	auto &pattern_here_grad = innerv?gp.pattern_grad_inner:gp.pattern_grad;
	auto &work_here_grad = innerv?gp.work_inner_grad:gp.work_grad;

	auto &pattern_here_hess = innerv?gp.pattern_hessian_inner:gp.pattern_hessian;
	auto &work_here_hess = innerv?gp.work_inner_hess:gp.work_hess;

	const auto &map_here = innerv?ctx->hessian_inner:ctx->hessian_outer;

	gp.fun.sparse_jac_rev(
		gp.x,
		pattern_here_grad,
		gp.unused_pattern,
		"cppad",
		work_here_grad);

	gp.fun.sparse_hes(
		gp.x,  
		gp.w,
		pattern_here_hess,              
		gp.unused_pattern,
		"cppad.symmetric",
		work_here_hess);

	const size_t NoutGrad = pattern_here_grad.nnz();
	const auto& cols = pattern_here_grad.col();
	const auto& vals = pattern_here_grad.val();

	for(size_t D=0; D < NoutGrad; D++) {
		out_grad[cols[D]] += vals[D];
	}

	const std::size_t start = map_here.map_p[ist];
	const std::size_t end   = map_here.map_p[ist+1];;
	const auto& vals_hess = pattern_here_hess.val();


	for(size_t Di=start;Di < end; Di++) {
		out_hess[ map_here.map_global[Di] ] += vals_hess[ map_here.map_local[Di] ];
	}
	return 0;
}





#endif
