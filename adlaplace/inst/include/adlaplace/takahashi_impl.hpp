#ifndef ADLAPLACE_TAKAHASHI_IMPL_HPP
#define ADLAPLACE_TAKAHASHI_IMPL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "adlaplace/chol_update_impl.hpp"

namespace adlaplace {
namespace chol {

// Forward-mode dual number so chol_update_csc / takahashi_selected_inv can
// deliver directional derivatives (used for second-order atomic passes).
struct Dual {
	double val;
	double dot;
	Dual() : val(0.0), dot(0.0) {}
	Dual(double v) : val(v), dot(0.0) {}
	Dual(double v, double d) : val(v), dot(d) {}
};

inline Dual operator+(const Dual& a, const Dual& b) {
	return Dual(a.val + b.val, a.dot + b.dot);
}
inline Dual operator-(const Dual& a, const Dual& b) {
	return Dual(a.val - b.val, a.dot - b.dot);
}
inline Dual operator-(const Dual& a) { return Dual(-a.val, -a.dot); }
inline Dual operator*(const Dual& a, const Dual& b) {
	return Dual(a.val * b.val, a.dot * b.val + a.val * b.dot);
}
inline Dual operator/(const Dual& a, const Dual& b) {
	const double v = a.val / b.val;
	return Dual(v, (a.dot - v * b.dot) / b.val);
}
inline Dual& operator+=(Dual& a, const Dual& b) {
	a.val += b.val;
	a.dot += b.dot;
	return a;
}
inline Dual& operator-=(Dual& a, const Dual& b) {
	a.val -= b.val;
	a.dot -= b.dot;
	return a;
}
inline Dual log(const Dual& a) { return Dual(std::log(a.val), a.dot / a.val); }

template <>
struct CholLogDetTraits<Dual> {
	static bool invalid_diag(const Dual& d) {
		return d.val <= 0.0 || !std::isfinite(d.val);
	}
	static Dual invalid_result() {
		const double nan = std::numeric_limits<double>::quiet_NaN();
		return Dual(nan, nan);
	}
	static Dual log_of(const Dual& d) { return log(d); }
};

// Look up a symmetric entry stored on a lower-triangular CSC pattern
// (column = min index, row = max index). Returns 0 when absent.
template <typename Scalar>
Scalar csc_sym_lower_entry(
	std::size_t r,
	std::size_t c,
	const std::vector<int>& p,
	const std::vector<int>& i,
	const std::vector<Scalar>& x)
{
	const std::size_t col = std::min(r, c);
	const std::size_t row = std::max(r, c);
	for (int pos = p[col]; pos < p[col + 1]; ++pos) {
		if (static_cast<std::size_t>(i[static_cast<std::size_t>(pos)]) == row) {
			return x[static_cast<std::size_t>(pos)];
		}
	}
	return Scalar(0);
}

// Takahashi selected inverse: given the LDL factor of P H P' from
// chol_update_csc (unit-lower L with pattern p_out/i_out and values x_out,
// diagonal d_out), compute Sigma = (P H P')^{-1} on the pattern of L
// (lower triangle including the diagonal). Recursion, for i >= j and
// columns processed from last to first:
//   Sigma_ij = delta_ij / d_j - sum_{k > j, L_kj != 0} L_kj * Sigma_ki
// The fill-in closure of the Cholesky pattern guarantees every Sigma_ki
// needed on the right-hand side lies on the pattern.
template <typename Scalar>
void takahashi_selected_inv(
	const std::vector<int>& p_out,
	const std::vector<int>& i_out,
	const std::vector<Scalar>& x_out,
	const std::vector<Scalar>& d_out,
	std::vector<Scalar>& sigma)
{
	const std::size_t n = p_out.size() - 1;
	sigma.assign(i_out.size(), Scalar(0));

	for (std::size_t jj = n; jj-- > 0;) {
		const std::size_t j = jj;
		int diag_pos = -1;
		// Off-diagonal entries of column j (any order: all Sigma_ki needed
		// live in columns > j, already computed).
		for (int pos = p_out[j]; pos < p_out[j + 1]; ++pos) {
			const std::size_t i =
				static_cast<std::size_t>(i_out[static_cast<std::size_t>(pos)]);
			if (i == j) {
				diag_pos = pos;
				continue;
			}
			Scalar acc = Scalar(0);
			for (int kpos = p_out[j]; kpos < p_out[j + 1]; ++kpos) {
				const std::size_t k = static_cast<std::size_t>(
					i_out[static_cast<std::size_t>(kpos)]);
				if (k == j) {
					continue; // unit diagonal of L
				}
				acc += x_out[static_cast<std::size_t>(kpos)] *
					csc_sym_lower_entry(k, i, p_out, i_out, sigma);
			}
			sigma[static_cast<std::size_t>(pos)] = -acc;
		}
		// Diagonal entry uses the off-diagonals of column j just computed
		// (stored at the same positions as L's column j).
		Scalar acc = Scalar(0);
		for (int kpos = p_out[j]; kpos < p_out[j + 1]; ++kpos) {
			const std::size_t k = static_cast<std::size_t>(
				i_out[static_cast<std::size_t>(kpos)]);
			if (k == j) {
				continue;
			}
			acc += x_out[static_cast<std::size_t>(kpos)] *
				sigma[static_cast<std::size_t>(kpos)];
		}
		if (diag_pos < 0) {
			Rcpp::stop("takahashi_selected_inv: missing diagonal in column %d",
			           static_cast<int>(j));
		}
		sigma[static_cast<std::size_t>(diag_pos)] = Scalar(1) / d_out[j] - acc;
	}
}

// Scatter the selected inverse (permuted ordering, pattern of L) onto the
// original-order CSC pattern of H: S_x[pos] = (H^{-1})_{Q_i[pos], col}.
// perm_inv maps original indices to permuted indices.
template <typename Scalar>
void selected_inv_scatter(
	const std::vector<int>& H_p,
	const std::vector<int>& H_i,
	const std::vector<int>& perm_inv,
	const std::vector<int>& p_out,
	const std::vector<int>& i_out,
	const std::vector<Scalar>& sigma,
	std::vector<Scalar>& S_x)
{
	const std::size_t n = H_p.size() - 1;
	S_x.assign(H_i.size(), Scalar(0));
	for (std::size_t c = 0; c < n; ++c) {
		for (int pos = H_p[c]; pos < H_p[c + 1]; ++pos) {
			const std::size_t r = static_cast<std::size_t>(
				H_i[static_cast<std::size_t>(pos)]);
			const std::size_t pr = static_cast<std::size_t>(perm_inv[r]);
			const std::size_t pc = static_cast<std::size_t>(perm_inv[c]);
			S_x[static_cast<std::size_t>(pos)] =
				csc_sym_lower_entry(pr, pc, p_out, i_out, sigma);
		}
	}
}

} // namespace chol
} // namespace adlaplace

#endif
