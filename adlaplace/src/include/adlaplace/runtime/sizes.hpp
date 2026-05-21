#ifndef API_SIZES_HPP
#define API_SIZES_HPP

#include <cppad/cppad.hpp>
#include <Rinternals.h>

#include <vector>
#include <cstddef>
#include <cstdio>

#include "adlaplace/api/backend.hpp"


// Get sparsity sizes for a single GroupPack
static int get_sizes(
    GroupPack& ad_pack, 
	int *n_inner,
	int *n_outer,
    int* nnz_grad_inner,
    int* nnz_grad_outer,
    int* nnz_hes_inner,
    int* nnz_hes_outer) {

	*n_inner = static_cast<int>(ad_pack.pattern_grad_inner.nc());	
	*n_outer = static_cast<int>(ad_pack.pattern_grad.nc());	

    *nnz_grad_outer = static_cast<int>(ad_pack.pattern_grad.nnz());
    *nnz_grad_inner = static_cast<int>(ad_pack.pattern_grad_inner.nnz());
    *nnz_hes_outer = static_cast<int>(ad_pack.pattern_hessian.nnz());
    *nnz_hes_inner = static_cast<int>(ad_pack.pattern_hessian_inner.nnz());

    return 0;
}

// Get sparsity pattern for a single GroupPack
static int get_pattern(
    GroupPack &ad_pack,
    int* pattern_grad_inner,
    int* pattern_grad_outer,
    int* pattern_hes_inner_row,
    int* pattern_hes_inner_col,
    int* pattern_hes_outer_row,
    int* pattern_hes_outer_col) {

    const size_t nnz_grad_inner = ad_pack.pattern_grad_inner.nnz();
    const auto& cols_grad_inner = ad_pack.pattern_grad_inner.col();
    for (size_t D = 0; D < nnz_grad_inner; D++) {
        pattern_grad_inner[D] = static_cast<int>(cols_grad_inner[D]);
    }
    const size_t nnz_grad_outer = ad_pack.pattern_grad.nnz();
    const auto& cols_grad_outer = ad_pack.pattern_grad.col();
    for (size_t D = 0; D < nnz_grad_outer; D++) {
        pattern_grad_outer[D] = static_cast<int>(cols_grad_outer[D]);
    }

    const size_t nnz_inner = ad_pack.pattern_hessian_inner.nnz();
    const auto& rows_hes_inner = ad_pack.pattern_hessian_inner.row();
    const auto& cols_hes_inner = ad_pack.pattern_hessian_inner.col();
    for (size_t D = 0; D < nnz_inner; D++) {
        pattern_hes_inner_row[D] = static_cast<int>(rows_hes_inner[D]);
        pattern_hes_inner_col[D] = static_cast<int>(cols_hes_inner[D]);
    }
    const size_t nnz_outer = ad_pack.pattern_hessian.nnz();
    const auto& rows_hes_outer = ad_pack.pattern_hessian.row();
    const auto& cols_hes_outer = ad_pack.pattern_hessian.col();
    for (size_t D = 0; D < nnz_outer; D++) {
        pattern_hes_outer_row[D] = static_cast<int>(rows_hes_outer[D]);
        pattern_hes_outer_col[D] = static_cast<int>(cols_hes_outer[D]);
    }

    return 0;
}

#endif
