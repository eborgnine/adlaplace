#include "adlaplace/chol_update.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include <Rcpp.h>
#include <Rinternals.h>
#include "adlaplace/backend.hpp"

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

void linv_update(
	const std::vector<int>& L_p,
	const std::vector<int>& L_i,
	const std::vector<double>& L_x,
	const std::vector<int>& Linv_p,
	const std::vector<int>& Linv_i,
	std::vector<double>& Linv_x)
{
	if (L_p.size() < 1) {
		Rcpp::stop("linv_update: L_p is empty");
	}
	const std::size_t n = L_p.size() - 1;
	if (Linv_p.size() != n + 1) {
		Rcpp::stop("linv_update: Linv_p length must be nrow + 1");
	}
	if (L_x.size() != L_i.size()) {
		Rcpp::stop("linv_update: L_i and L_x length must match");
	}
	Linv_x.assign(Linv_i.size(), 0.0);

	std::vector<double> linv_col(n, 0.0);

	for (std::size_t j = 0; j < n; ++j) {
		std::fill(linv_col.begin(), linv_col.end(), 0.0);
		for (std::size_t pos = static_cast<std::size_t>(Linv_p[j]);
		     pos < static_cast<std::size_t>(Linv_p[j + 1]);
		     ++pos) {
			const std::size_t i = static_cast<std::size_t>(Linv_i[pos]);
			if (i == j) {
				linv_col[i] = 1.0;
				Linv_x[pos] = 1.0;
			}
		}

		for (std::size_t pos = static_cast<std::size_t>(Linv_p[j]);
		     pos < static_cast<std::size_t>(Linv_p[j + 1]);
		     ++pos) {
			const std::size_t i = static_cast<std::size_t>(Linv_i[pos]);
			if (i <= j) {
				continue;
			}

			double sum = 0.0;
			for (std::size_t k = j; k < i; ++k) {
				double l_ik = 0.0;
				for (std::size_t lpos = static_cast<std::size_t>(L_p[k]);
				     lpos < static_cast<std::size_t>(L_p[k + 1]);
				     ++lpos) {
					if (static_cast<std::size_t>(L_i[lpos]) == i) {
						l_ik = L_x[lpos];
						break;
					}
				}
				if (l_ik != 0.0) {
					sum += l_ik * linv_col[k];
				}
			}
			Linv_x[pos] = -sum;
			linv_col[i] = Linv_x[pos];
		}
	}
}

double csc_entry(
	const std::size_t row,
	const std::size_t col,
	const std::vector<int>& p,
	const std::vector<int>& i,
	const std::vector<double>& x)
{
	for (std::size_t pos = static_cast<std::size_t>(p[col]);
	     pos < static_cast<std::size_t>(p[col + 1]);
	     ++pos) {
		if (static_cast<std::size_t>(i[pos]) == row) {
			return x[pos];
		}
	}
	return 0.0;
}

void half_h_inv_update(
	const std::vector<int>& Linv_p,
	const std::vector<int>& Linv_i,
	const std::vector<double>& Linv_x,
	const std::vector<double>& d,
	const std::vector<int>& perm_inv,
	const std::vector<int>& half_H_inv_p,
	const std::vector<int>& half_H_inv_i,
	std::vector<double>& half_H_inv_x)
{
	if (Linv_p.size() < 1) {
		Rcpp::stop("half_h_inv_update: Linv_p is empty");
	}
	const std::size_t n = Linv_p.size() - 1;
	if (Linv_x.size() != Linv_i.size()) {
		Rcpp::stop("half_h_inv_update: Linv_i and Linv_x length must match");
	}
	if (d.size() != n) {
		Rcpp::stop("half_h_inv_update: d length must equal nrow");
	}
	if (perm_inv.size() != n) {
		Rcpp::stop("half_h_inv_update: perm_inv length must equal nrow");
	}
	if (half_H_inv_p.size() != n + 1) {
		Rcpp::stop("half_h_inv_update: half_H_inv_p length must be nrow + 1");
	}
	if (half_H_inv_x.size() != half_H_inv_i.size()) {
		half_H_inv_x.assign(half_H_inv_i.size(), 0.0);
	} else {
		std::fill(half_H_inv_x.begin(), half_H_inv_x.end(), 0.0);
	}

	for (std::size_t col = 0; col < n; ++col) {
		for (std::size_t pos = static_cast<std::size_t>(half_H_inv_p[col]);
		     pos < static_cast<std::size_t>(half_H_inv_p[col + 1]);
		     ++pos) {
			const std::size_t row = static_cast<std::size_t>(half_H_inv_i[pos]);
			const std::size_t perm_row = static_cast<std::size_t>(perm_inv[row]);
			const double scale = std::pow(d[col], -0.5);
			half_H_inv_x[pos] = csc_entry(col, perm_row, Linv_p, Linv_i, Linv_x) * scale;
		}
	}
}

void h_inv_update(
	const std::vector<int>& half_H_inv_p,
	const std::vector<int>& half_H_inv_i,
	const std::vector<double>& half_H_inv_x,
	const std::vector<int>& H_inv_p,
	const std::vector<int>& H_inv_i,
	std::vector<double>& H_inv_x)
{
	if (half_H_inv_p.size() < 1) {
		Rcpp::stop("h_inv_update: half_H_inv_p is empty");
	}
	const std::size_t n = half_H_inv_p.size() - 1;
	if (half_H_inv_x.size() != half_H_inv_i.size()) {
		Rcpp::stop("h_inv_update: half_H_inv_i and half_H_inv_x length must match");
	}
	if (H_inv_p.size() != n + 1) {
		Rcpp::stop("h_inv_update: H_inv_p length must be nrow + 1");
	}
	if (H_inv_x.size() != H_inv_i.size()) {
		H_inv_x.assign(H_inv_i.size(), 0.0);
	} else {
		std::fill(H_inv_x.begin(), H_inv_x.end(), 0.0);
	}

	std::vector<std::vector<std::pair<int, double>>> rows_by_col(n);
	for (std::size_t col = 0; col < n; ++col) {
		for (std::size_t pos = static_cast<std::size_t>(half_H_inv_p[col]);
		     pos < static_cast<std::size_t>(half_H_inv_p[col + 1]);
		     ++pos) {
			rows_by_col[col].emplace_back(
				half_H_inv_i[pos],
				half_H_inv_x[pos]
			);
		}
	}

	std::vector<std::vector<std::pair<int, double>>> row_entries(n);
	for (std::size_t col = 0; col < n; ++col) {
		for (const auto& entry : rows_by_col[col]) {
			row_entries[static_cast<std::size_t>(entry.first)].emplace_back(
				static_cast<int>(col),
				entry.second
			);
		}
	}

	for (std::size_t col = 0; col < n; ++col) {
		for (std::size_t pos = static_cast<std::size_t>(H_inv_p[col]);
		     pos < static_cast<std::size_t>(H_inv_p[col + 1]);
		     ++pos) {
			const std::size_t row = static_cast<std::size_t>(H_inv_i[pos]);
			if (row > col) {
				continue;
			}
			const auto& a = row_entries[row];
			const auto& b = row_entries[col];
			std::size_t ia = 0;
			std::size_t ib = 0;
			double sum = 0.0;
			while (ia < a.size() && ib < b.size()) {
				if (a[ia].first < b[ib].first) {
					++ia;
				} else if (a[ia].first > b[ib].first) {
					++ib;
				} else {
					sum += a[ia].second * b[ib].second;
					++ia;
					++ib;
				}
			}
			H_inv_x[pos] = sum;
		}
	}
}

namespace {

void copy_csc_pi_from_s4(
	const Rcpp::S4& mat,
	std::vector<int>& p_out,
	std::vector<int>& i_out)
{
	const Rcpp::IntegerVector p = mat.slot("p");
	const Rcpp::IntegerVector i = mat.slot("i");
	p_out.assign(p.begin(), p.end());
	i_out.assign(i.begin(), i.end());
}

CholPattern chol_pattern_from_list(const Rcpp::List& cil, const int n_gamma) {
	CholPattern pattern;
	pattern.n = n_gamma;
	if (n_gamma == 0) {
		return pattern;
	}

	if (!cil.containsElementNamed("L1") ||
	    !cil.containsElementNamed("Linv") ||
	    !cil.containsElementNamed("perm") ||
	    !cil.containsElementNamed("perm_inv") ||
	    !cil.containsElementNamed("half_H_inv") ||
	    !cil.containsElementNamed("H_inv")) {
		Rcpp::stop(
			"chol_inner_list must contain L1, Linv, perm, perm_inv, half_H_inv, and H_inv"
		);
	}

	const Rcpp::S4 L1 = cil["L1"];
	const Rcpp::S4 Linv = cil["Linv"];
	const Rcpp::S4 half_H_inv = cil["half_H_inv"];
	const Rcpp::S4 H_inv = cil["H_inv"];
	const Rcpp::IntegerVector perm_r = cil["perm"];
	const Rcpp::IntegerVector perm_inv_r = cil["perm_inv"];

	const Rcpp::IntegerVector L1_dim = L1.slot("Dim");
	if (L1_dim.size() < 1 || static_cast<int>(L1_dim[0]) != n_gamma) {
		Rcpp::stop(
			"chol_inner_list$L1 nrow (%d) must equal Ngamma (%d)",
			static_cast<int>(L1_dim[0]),
			n_gamma
		);
	}
	if (perm_r.size() != n_gamma || perm_inv_r.size() != n_gamma) {
		Rcpp::stop("chol_inner_list perm vectors must have length Ngamma");
	}

	copy_csc_pi_from_s4(L1, pattern.L1_p, pattern.L1_i);
	copy_csc_pi_from_s4(Linv, pattern.Linv_p, pattern.Linv_i);
	copy_csc_pi_from_s4(half_H_inv, pattern.half_H_inv_p, pattern.half_H_inv_i);
	copy_csc_pi_from_s4(H_inv, pattern.H_inv_p, pattern.H_inv_i);
	if (cil.containsElementNamed("trace_columns")) {
		const Rcpp::S4 trace_columns = cil["trace_columns"];
		copy_csc_pi_from_s4(
			trace_columns,
			pattern.trace_columns_p,
			pattern.trace_columns_i
		);
	}

	pattern.perm.resize(static_cast<std::size_t>(n_gamma));
	pattern.perm_inv.resize(static_cast<std::size_t>(n_gamma));
	for (int j = 0; j < n_gamma; ++j) {
		// hessian_map() stores 0-based perm / perm_inv (inner Hessian uses index1 = FALSE).
		pattern.perm[static_cast<std::size_t>(j)] = perm_r[j];
		pattern.perm_inv[static_cast<std::size_t>(j)] = perm_inv_r[j];
	}

	return pattern;
}

} // namespace

void ad_fun_attach_chol_pattern_from_list(ad_fun& shards, const Rcpp::List& hessian_pack) {
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

	if (!hessian_pack.containsElementNamed("chol_inner_list")) {
		Rcpp::stop("hessian_pack must contain chol_inner_list from hessian_map()");
	}

	const Rcpp::List cil = hessian_pack["chol_inner_list"];
	if (cil.size() == 0) {
		Rcpp::stop(
			"chol_inner_list is empty but Ngamma = %d; symbolic LDL pattern is required",
			inner_n
		);
	}

	shards.chol_pattern = chol_pattern_from_list(cil, inner_n);
	if (shards.chol_pattern.n != inner_n) {
		Rcpp::stop(
			"inner chol pattern dimension (%d) does not match inner Hessian template (%d)",
			shards.chol_pattern.n,
			inner_n
		);
	}
}
