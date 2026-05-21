#include "chol_update.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <Rcpp.h>
#include <Rinternals.h>

#include "adlaplace/api/backend.hpp"

namespace {

// Upper-triangle CSC (Eigen column-major): entry (row, col) in column max(row, col).
double sym_entry(
	size_t row,
	size_t col,
	const Eigen::SparseMatrix<double>& A)
{
	const size_t c = std::max(row, col);
	const size_t r = std::min(row, col);
	for (Eigen::SparseMatrix<double>::InnerIterator it(A, static_cast<int>(c)); it; ++it) {
		if (static_cast<size_t>(it.row()) == r) {
			return it.value();
		}
	}
	return 0.0;
}

// (H_perm)_{r,j} = H_{perm[r], perm[j]} with Matrix perm[j] = original index for permuted j.
double sym_entry_perm(
	size_t row,
	size_t col,
	const Eigen::SparseMatrix<double>& H,
	const std::vector<int>& perm)
{
	return sym_entry(
		static_cast<size_t>(perm[row]),
		static_cast<size_t>(perm[col]),
		H
	);
}

std::vector<int> copy_int_slot(const Rcpp::S4& obj, const char* name) {
	return Rcpp::as<std::vector<int>>(obj.slot(name));
}

} // namespace

// Numeric LDL of P H P' in permuted ordering; H is original (upper CSC), perm from chol_inner.
// Fills unit-lower L (off-diagonal in x_out, pattern p_out/i_out) and diagonal D.
// Returns log(det) = sum(log(D)); NA_REAL if factorization D is invalid.
double chol_update(
	const Eigen::SparseMatrix<double>& H,
	const std::vector<int>& perm,
	const std::vector<int>& p_out,
	const std::vector<int>& i_out,
	std::vector<double>& x_out,
	std::vector<double>& d_out)
{
	if (H.rows() != H.cols()) {
		Rcpp::stop("chol_update: H must be square");
	}
	const size_t n = static_cast<size_t>(H.rows());
	if (perm.size() != n) {
		Rcpp::stop(
			"chol_update: perm length (%d) must match H nrow (%d)",
			static_cast<int>(perm.size()),
			static_cast<int>(n)
		);
	}
	if (p_out.size() != n + 1) {
		Rcpp::stop("chol_update: p_out length must be nrow + 1");
	}
	if (d_out.size() != n) {
		d_out.assign(n, 0.0);
	}
	if (x_out.size() != i_out.size()) {
		x_out.assign(i_out.size(), 0.0);
	}

	std::vector<double> y(n, 0.0);

	for (size_t j = 0; j < n; ++j) {
		for (size_t r = 0; r < n; ++r) {
			y[r] = sym_entry_perm(r, j, H, perm);
		}

		for (size_t k = 0; k < j; ++k) {
			bool has_l_jk = false;
			double l_jk = 0.0;
			for (size_t pos = p_out[k]; pos < p_out[k + 1]; ++pos) {
				if (static_cast<size_t>(i_out[pos]) == j) {
					l_jk = x_out[pos];
					has_l_jk = true;
					break;
				}
			}
			if (!has_l_jk) {
				continue;
			}
			const double alpha = l_jk * d_out[k];
			for (size_t pos = p_out[k]; pos < p_out[k + 1]; ++pos) {
				const size_t row = static_cast<size_t>(i_out[pos]);
				if (row < j) {
					continue;
				}
				y[row] -= alpha * x_out[pos];
			}
		}

		d_out[j] = y[j];
		for (size_t pos = p_out[j]; pos < p_out[j + 1]; ++pos) {
			const size_t row = static_cast<size_t>(i_out[pos]);
			if (row == j) {
				x_out[pos] = 1.0;
			} else if (row > j) {
				x_out[pos] = y[row] / d_out[j];
			}
		}
	}

	double log_det = 0.0;
	for (size_t j = 0; j < n; ++j) {
		const double dj = d_out[j];
		if (dj <= 0.0 || !std::isfinite(dj)) {
			return NA_REAL;
		}
		log_det += std::log(dj);
	}
	return log_det;
}

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
	pattern.perm = copy_int_slot(obj, "perm");

	const Rcpp::IntegerVector dim = obj.slot("Dim");
	if (dim.size() != 2) {
		Rcpp::stop("chol_inner Dim slot must have length 2");
	}
	const int n = dim[0];
	const int ncol = dim[1];
	if (n != ncol || n < 0) {
		Rcpp::stop("chol_inner must be a square matrix (Dim[0] == Dim[1] >= 0)");
	}
	if (pattern.perm.size() != static_cast<size_t>(n)) {
		Rcpp::stop(
			"chol_inner perm length (%d) does not match dimension (%d)",
			static_cast<int>(pattern.perm.size()),
			n
		);
	}
	for (size_t j = 0; j < static_cast<size_t>(n); ++j) {
		const int orig = pattern.perm[j];
		if (orig < 0 || orig >= n) {
			Rcpp::stop(
				"chol_inner perm[%d] = %d out of range [0, %d]",
				static_cast<int>(j),
				orig,
				n - 1
			);
		}
	}

	pattern.n = n;
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
