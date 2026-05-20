#include "chol_update.hpp"

#include <Rinternals.h>

#include "adlaplace/runtime/backend.hpp"

namespace {

std::vector<int> copy_int_slot(const Rcpp::S4& obj, const char* name) {
	return Rcpp::as<std::vector<int>>(obj.slot(name));
}

// Matrix perm[j] = original row for permuted row j; Eigen wants inverse in pinv.
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pinv_from_perm(
	const std::vector<int>& perm,
	int n)
{
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> pinv(n);
	for (int j = 0; j < n; ++j) {
		const int orig = perm[static_cast<size_t>(j)];
		if (orig < 0 || orig >= n) {
			Rcpp::stop(
				"chol_pattern perm[%d] = %d out of range [0, %d]",
				j,
				orig,
				n - 1
			);
		}
		pinv.indices()[orig] = j;
	}
	return pinv;
}

} // namespace

CholPattern chol_pattern_from_sexp(SEXP chol_inner) {
	if (Rf_isNull(chol_inner)) {
		Rcpp::stop("chol_inner must not be NULL");
	}
	if (!Rf_inherits(chol_inner, "dCHMsimpl")) {
		Rcpp::stop("chol_inner must be a dCHMsimpl object");
	}

	CholPattern pattern;
	const Rcpp::S4 obj(chol_inner);
	pattern.p = copy_int_slot(obj, "p");
	pattern.i = copy_int_slot(obj, "i");
	pattern.nz = copy_int_slot(obj, "nz");
	pattern.nxt = copy_int_slot(obj, "nxt");
	pattern.prv = copy_int_slot(obj, "prv");
	pattern.type = copy_int_slot(obj, "type");
	pattern.colcount = copy_int_slot(obj, "colcount");
	pattern.perm = copy_int_slot(obj, "perm");
	pattern.dim = copy_int_slot(obj, "Dim");

	if (pattern.dim.size() != 2) {
		Rcpp::stop("chol_inner Dim slot must have length 2");
	}
	const int n = pattern.dim[0];
	const int ncol = pattern.dim[1];
	if (n != ncol || n < 0) {
		Rcpp::stop("chol_inner must be a square matrix (Dim[0] == Dim[1] >= 0)");
	}
	if (static_cast<int>(pattern.perm.size()) != n) {
		Rcpp::stop(
			"chol_inner perm length (%d) does not match dimension (%d)",
			static_cast<int>(pattern.perm.size()),
			n
		);
	}

	pattern.n = n;
	pattern.pinv = pinv_from_perm(pattern.perm, n);
	return pattern;
}

void ad_groups_attach_chol_pattern(ad_groups& groups, const Rcpp::List& ad_fun) {
	if (groups.chol_pattern.n > 0) {
		return;
	}
	if (!groups.hessians_attached) {
		Rcpp::stop("attach Hessian templates before chol_inner");
	}

	if (!ad_fun.containsElementNamed("chol_inner") || Rf_isNull(ad_fun["chol_inner"])) {
		Rcpp::stop("ad_fun must contain chol_inner from hessian_map()");
	}

	const int inner_n = static_cast<int>(groups.hessian_inner.rows());
	groups.chol_pattern = chol_pattern_from_sexp(ad_fun["chol_inner"]);
	if (groups.chol_pattern.n != inner_n) {
		Rcpp::stop(
			"chol_inner dimension (%d) does not match inner Hessian template (%d)",
			groups.chol_pattern.n,
			inner_n
		);
	}
}
