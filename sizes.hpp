#ifndef API_SIZES_HPP
#define API_SIZES_HPP

#include <cppad/cppad.hpp>
#include <Rinternals.h>

#include <vector>
#include <cstddef>
#include <cstdio>

#include "adlaplace/runtime/backend.hpp"


static int get_sizes(
    GroupPack& ad_pack, 
    int* Nparams, 
	int* Nbeta, int* Ngamma, int* Ntheta,
	int* nnz_grad_inner,
	int* nnz_grad_outer,
	int* nnz_hes_inner,
	int* nnz_hes_outer) {

    *nnz_grad_outer = ad_pack.pattern_grad_outer.nnz();
	*nnz_grad_inner = ad_pack.pattern_grad_inner.nnz();
	*nnz_hes_outer = ad_pack.pattern_hes_outer.nnz();
	*nnz_hes_inner = ad_pack.pattern_hes_inner.nnz();
    
    *Nbeta = ad_pack.Nbeta;
    *Ngamma = ad_pack.Ngamma;
    *Ntheta= ad_pack.Ntheta;
    *Nparams = ad_pack.x.size();

	return 0;
}

static int get_pattern(
	GroupPack &ad_pack,
	int* pattern_grad_inner,
	int* pattern_grad_outer,
	int* pattern_hes_inner_row,
	int* pattern_hes_inner_col,
	int* pattern_hes_outer_row,
	int* pattern_hes_outer_col) {

	const size_t nnz_grad_inner = ad_pack.pattern_grad_inner.nnz();
	for (size_t D = 0; D < nnz_grad_inner; D++) {
		pattern_grad_inner[D] = ad_pack.pattern_grad_inner.col[D];
	}
	const size_t nnz_grad_outer = ad_pack.pattern_grad_outer.nnz();
	for (size_t D = 0; D < nnz_grad_outer; D++) {
		pattern_grad_outer[D] = ad_pack.pattern_grad_outer.col[D];
	}

	const size_t nnz_inner = ad_pack.pattern_hes_inner.nnz();
	for (size_t D = 0; D < nnz_inner; D++) {
		pattern_hes_inner_row[D] = ad_pack.pattern_hes_inner.row[D];
		pattern_hes_inner_col[D] = ad_pack.pattern_hes_inner.col[D];
	}
	const size_t nnz_outer = ad_pack.pattern_hes.nnz();
	for (size_t D = 0; D < nnz_outer; D++) {
		pattern_hes_outer_row[D] = ad_pack.pattern_hes_outer.row[D];
		pattern_hes_outer_col[D] = ad_pack.pattern_hes_outer.col[D];
	}
}


#endif