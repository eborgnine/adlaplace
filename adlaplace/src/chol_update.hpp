#ifndef ADLAPLACE_CHOL_UPDATE_HPP
#define ADLAPLACE_CHOL_UPDATE_HPP

#include <Rcpp.h>
#include <Eigen/Sparse>

#include "adlaplace/api/backend.hpp"

#include <vector>

// Returns log(det(H_perm)) = sum(log(D)); NA_REAL if any D is non-positive or non-finite.
double chol_update(
	const Eigen::SparseMatrix<double>& H,
	const std::vector<int>& perm,
	const std::vector<int>& p_out,
	const std::vector<int>& i_out,
	std::vector<double>& x_out,
	std::vector<double>& d_out
);

CholPattern chol_pattern_from_sexp(SEXP chol_inner);

CholPattern chol_pattern_from_list(const Rcpp::List& chol_inner_list);

void ad_groups_attach_chol_pattern(ad_groups& groups, const Rcpp::List& ad_fun);

#endif
