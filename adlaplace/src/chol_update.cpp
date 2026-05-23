#include "chol_update.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <Rcpp.h>
#include <Rinternals.h>
#include <Eigen/SparseCholesky>

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

Eigen::SparseMatrix<double> dummy_spd_from_template(const hessian_template& inner) {
	const int n = static_cast<int>(inner.rows());
	Eigen::SparseMatrix<double> H(n, n);
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(static_cast<size_t>(inner.nonZeros()) * 2);
	for (int col = 0; col < inner.outerSize(); ++col) {
		for (hessian_template::InnerIterator it(inner, col); it; ++it) {
			const int row = static_cast<int>(it.row());
			if (row > col) {
				continue;
			}
			const double v = (row == col) ? 10.0 : 1.0;
			triplets.emplace_back(row, col, v);
			if (row != col) {
				triplets.emplace_back(col, row, v);
			}
		}
	}
	H.setFromTriplets(triplets.begin(), triplets.end());
	return H;
}

void copy_lower_csc_pattern(
	const Eigen::SparseMatrix<double, Eigen::ColMajor>& L,
	std::vector<int>& p_out,
	std::vector<int>& i_out)
{
	Eigen::SparseMatrix<double, Eigen::ColMajor> Lc = L;
	Lc.makeCompressed();
	const int n = static_cast<int>(Lc.cols());
	p_out.assign(static_cast<size_t>(n + 1), 0);
	i_out.clear();
	i_out.reserve(static_cast<size_t>(Lc.nonZeros()));
	for (int col = 0; col < n; ++col) {
		p_out[static_cast<size_t>(col)] = static_cast<int>(i_out.size());
		for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(Lc, col); it; ++it) {
			const int row = static_cast<int>(it.row());
			if (row >= col) {
				i_out.push_back(row);
			}
		}
	}
	p_out[static_cast<size_t>(n)] = static_cast<int>(i_out.size());
}

CholPattern chol_pattern_identity_lower_from_template(const hessian_template& inner) {
	CholPattern pattern;
	const int n = static_cast<int>(inner.rows());
	pattern.n = n;
	if (n == 0) {
		return pattern;
	}
	pattern.perm.assign(static_cast<size_t>(n), 0);
	pattern.perm_inv.assign(static_cast<size_t>(n), 0);
	for (int j = 0; j < n; ++j) {
		pattern.perm[static_cast<size_t>(j)] = j;
		pattern.perm_inv[static_cast<size_t>(j)] = j;
	}
	pattern.L1_p.assign(static_cast<size_t>(n + 1), 0);
	for (int col = 0; col < n; ++col) {
		for (hessian_template::InnerIterator it(inner, col); it; ++it) {
			const int row = static_cast<int>(it.row());
			if (row >= col) {
				pattern.L1_i.push_back(row);
			}
		}
		pattern.L1_p[static_cast<size_t>(col + 1)] =
			static_cast<int>(pattern.L1_i.size());
	}
	pattern.Linv_p = pattern.L1_p;
	pattern.Linv_i = pattern.L1_i;
	return pattern;
}

} // namespace

// Numeric LDL of P H P' in permuted ordering; H is original (upper CSC), perm 0-based.
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

CholPattern chol_pattern_from_inner_template(const hessian_template& inner) {
	const int n = static_cast<int>(inner.rows());
	if (n == 0) {
		return CholPattern();
	}

	const Eigen::SparseMatrix<double> H = dummy_spd_from_template(inner);
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(H);
	if (ldlt.info() != Eigen::Success) {
		return chol_pattern_identity_lower_from_template(inner);
	}

	CholPattern pattern;
	pattern.n = n;

	const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P = ldlt.permutationP();
	pattern.perm.resize(static_cast<size_t>(n));
	pattern.perm_inv.assign(static_cast<size_t>(n), 0);
	for (int j = 0; j < n; ++j) {
		const int orig = static_cast<int>(P.indices()[j]);
		pattern.perm[static_cast<size_t>(j)] = orig;
		pattern.perm_inv[static_cast<size_t>(orig)] = j;
	}

	const Eigen::SparseMatrix<double, Eigen::ColMajor> L = ldlt.matrixL();
	copy_lower_csc_pattern(L, pattern.L1_p, pattern.L1_i);
	pattern.Linv_p = pattern.L1_p;
	pattern.Linv_i = pattern.L1_i;
	return pattern;
}

void ad_fun_attach_chol_pattern_from_template(ad_fun& shards) {
	if (shards.chol_pattern.n > 0) {
		return;
	}
	if (!shards.hessians_attached) {
		Rcpp::stop("attach Hessian templates before chol pattern");
	}
	const int inner_n = static_cast<int>(shards.hessian_inner.rows());
	if (inner_n == 0) {
		shards.chol_pattern = CholPattern();
		return;
	}
	shards.chol_pattern = chol_pattern_from_inner_template(shards.hessian_inner);
	if (shards.chol_pattern.n != inner_n) {
		Rcpp::stop(
			"inner chol pattern dimension (%d) does not match inner Hessian template (%d)",
			shards.chol_pattern.n,
			inner_n
		);
	}
}
