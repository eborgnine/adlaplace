#ifndef ADLAPLACEFEM_TAKAHASHI_DAVIS_HPP
#define ADLAPLACEFEM_TAKAHASHI_DAVIS_HPP

// Takahashi–Davis selected inverse.
//
// Given the unit-lower CSC factor L and the LDL diagonal D of a symmetric
// positive-definite matrix (diagonal of L explicit and equal to one; row
// indices in each column sorted), write (L D L')^{-1} on the sparsity
// pattern of L. Columns run from last to first. Already computed entries
// are scattered into a work vector of length n and the Takahashi update
// reads that vector by row index.
//
// Header-only. A package that lists LinkingTo: adlaplaceFem compiles this
// template into its own shared library:
//
//   #include <adlaplaceFem/takahashi_davis.hpp>
//   adlaplaceFem::takahashi_davis(Lp, Li, Lx, D, sigma);
//
// A caller that has only Q factors it first (chol_update_csc in adlaplace).

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace adlaplaceFem {

template <typename Scalar>
void takahashi_davis(
    const std::vector<int>& Lp,
    const std::vector<int>& Li,
    const std::vector<Scalar>& Lx,
    const std::vector<Scalar>& D,
    std::vector<Scalar>& sigma) {
	if (Lp.size() < 2) {
		throw std::invalid_argument("takahashi_davis: L has no columns");
	}
	const std::size_t n = Lp.size() - 1;
	if (D.size() != n) {
		throw std::invalid_argument("takahashi_davis: D length must equal nrow(L)");
	}
	if (Lx.size() != Li.size()) {
		throw std::invalid_argument("takahashi_davis: Lx and Li lengths differ");
	}
	sigma.assign(Li.size(), Scalar(0));

	std::vector<int> diag_pos(n, -1);
	// Off-diagonal L(row, col) positions, grouped by row, so column j can
	// walk the k < j with L(j, k) != 0 without scanning earlier columns.
	struct Link {
		int col;
		int pos;
	};
	std::vector<std::vector<Link>> links(n);
	for (std::size_t k = 0; k < n; ++k) {
		for (int pos = Lp[k]; pos < Lp[k + 1]; ++pos) {
			const int row = Li[static_cast<std::size_t>(pos)];
			if (row == static_cast<int>(k)) {
				diag_pos[k] = pos;
				sigma[static_cast<std::size_t>(pos)] = Scalar(1) / D[k];
			} else {
				links[static_cast<std::size_t>(row)].push_back(
				    Link{static_cast<int>(k), pos});
			}
		}
		if (diag_pos[k] < 0) {
			throw std::invalid_argument(
			    "takahashi_davis: missing diagonal in L");
		}
	}
	for (std::size_t row = 0; row < n; ++row) {
		auto& hit = links[row];
		// Decreasing column index: z[i] for i > k is already known.
		std::sort(hit.begin(), hit.end(), [](const Link& a, const Link& b) {
			return a.col > b.col;
		});
	}

	std::vector<Scalar> z(n, Scalar(0));
	for (std::size_t jj = n; jj-- > 0;) {
		const std::size_t j = jj;
		for (int pos = Lp[j]; pos < Lp[j + 1]; ++pos) {
			z[static_cast<std::size_t>(Li[static_cast<std::size_t>(pos)])] =
			    sigma[static_cast<std::size_t>(pos)];
		}

		const std::vector<Link>& hit = links[j];
		for (std::size_t t = 0; t < hit.size(); ++t) {
			const int k = hit[t].col;
			Scalar zkj(0);
			for (int p = Lp[static_cast<std::size_t>(k)];
			     p < Lp[static_cast<std::size_t>(k) + 1]; ++p) {
				const int i = Li[static_cast<std::size_t>(p)];
				if (i > k) {
					zkj -= Lx[static_cast<std::size_t>(p)] *
					       z[static_cast<std::size_t>(i)];
				}
			}
			z[static_cast<std::size_t>(k)] = zkj;
		}

		for (std::size_t t = 0; t < hit.size(); ++t) {
			const int k = hit[t].col;
			const Scalar ljk = Lx[static_cast<std::size_t>(hit[t].pos)];
			const std::size_t ks = static_cast<std::size_t>(k);
			for (int p = Lp[ks]; p < Lp[ks + 1]; ++p) {
				const std::size_t i =
				    static_cast<std::size_t>(Li[static_cast<std::size_t>(p)]);
				sigma[static_cast<std::size_t>(p)] -= z[i] * ljk;
			}
		}

		for (std::size_t t = 0; t < hit.size(); ++t) {
			z[static_cast<std::size_t>(hit[t].col)] = Scalar(0);
		}
		for (int pos = Lp[j]; pos < Lp[j + 1]; ++pos) {
			const std::size_t i =
			    static_cast<std::size_t>(Li[static_cast<std::size_t>(pos)]);
			sigma[static_cast<std::size_t>(pos)] = z[i];
			z[i] = Scalar(0);
		}
	}
}

}  // namespace adlaplaceFem

#endif
