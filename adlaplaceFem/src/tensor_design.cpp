#include <Rcpp.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

// Serial tensor-product fill. Row r of the design is the Kronecker product of
// the nonzero B-spline weights at that cell: column i + (j - 1) * nx gets
// Bx[r, i] * By[r, j]. Triplets are 1-based, as in Matrix::summary().
// [[Rcpp::export]]
Rcpp::S4 tensor_design_triplets(Rcpp::IntegerVector bx_i, Rcpp::IntegerVector bx_j,
                                Rcpp::NumericVector bx_x, Rcpp::IntegerVector by_i,
                                Rcpp::IntegerVector by_j, Rcpp::NumericVector by_x,
                                int n, int nx, int ny) {
	if (bx_i.size() != bx_j.size() || bx_i.size() != bx_x.size()) {
		Rcpp::stop("tensor_design: Bx triplet lengths differ");
	}
	if (by_i.size() != by_j.size() || by_i.size() != by_x.size()) {
		Rcpp::stop("tensor_design: By triplet lengths differ");
	}
	if (n < 0 || nx < 0 || ny < 0) {
		Rcpp::stop("tensor_design: dimensions must be non-negative");
	}
	const long long ncol_ll =
	    static_cast<long long>(nx) * static_cast<long long>(ny);
	if (ncol_ll > static_cast<long long>(std::numeric_limits<int>::max())) {
		Rcpp::stop("tensor_design: basis dimension exceeds integer range");
	}
	const int ncol = static_cast<int>(ncol_ll);

	struct Weight {
		int j;
		double x;
	};
	std::vector<std::vector<Weight>> bx_rows(static_cast<std::size_t>(n));
	std::vector<std::vector<Weight>> by_rows(static_cast<std::size_t>(n));

	auto bucket = [&](const Rcpp::IntegerVector& rows, const Rcpp::IntegerVector& cols,
	                  const Rcpp::NumericVector& vals,
	                  std::vector<std::vector<Weight>>& dest, int n_basis,
	                  const char* which) {
		const R_xlen_t m = rows.size();
		for (R_xlen_t k = 0; k < m; ++k) {
			if (Rcpp::IntegerVector::is_na(rows[k]) || Rcpp::IntegerVector::is_na(cols[k])) {
				Rcpp::stop("tensor_design: %s triplet has a missing index", which);
			}
			const int r = rows[k];
			const int j = cols[k];
			if (r < 1 || r > n || j < 1 || j > n_basis) {
				Rcpp::stop("tensor_design: %s triplet index out of range", which);
			}
			dest[static_cast<std::size_t>(r - 1)].push_back(Weight{j, vals[k]});
		}
	};
	bucket(bx_i, bx_j, bx_x, bx_rows, nx, "Bx");
	bucket(by_i, by_j, by_x, by_rows, ny, "By");

	struct Trip {
		int row;
		int col;
		double val;
		bool operator<(const Trip& other) const {
			if (col != other.col) {
				return col < other.col;
			}
			return row < other.row;
		}
	};

	std::size_t nnz_bound = 0;
	for (int r = 0; r < n; ++r) {
		nnz_bound += bx_rows[static_cast<std::size_t>(r)].size() *
		             by_rows[static_cast<std::size_t>(r)].size();
	}
	std::vector<Trip> trips;
	trips.reserve(nnz_bound);
	for (int r = 0; r < n; ++r) {
		const std::vector<Weight>& bx = bx_rows[static_cast<std::size_t>(r)];
		const std::vector<Weight>& by = by_rows[static_cast<std::size_t>(r)];
		if (bx.empty() || by.empty()) {
			continue;
		}
		for (const Weight& a : bx) {
			for (const Weight& b : by) {
				const int col = (a.j - 1) + (b.j - 1) * nx;
				trips.push_back(Trip{r, col, a.x * b.x});
			}
		}
	}

	std::sort(trips.begin(), trips.end());
	std::vector<Trip> uniq;
	uniq.reserve(trips.size());
	for (const Trip& t : trips) {
		if (!uniq.empty() && uniq.back().col == t.col && uniq.back().row == t.row) {
			uniq.back().val += t.val;
		} else {
			uniq.push_back(t);
		}
	}
	if (uniq.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
		Rcpp::stop("tensor_design: too many nonzeros");
	}

	const int nnz = static_cast<int>(uniq.size());
	Rcpp::IntegerVector i(nnz);
	Rcpp::NumericVector x(nnz);
	Rcpp::IntegerVector p(ncol + 1);
	for (int k = 0; k < nnz; ++k) {
		i[k] = uniq[static_cast<std::size_t>(k)].row;
		x[k] = uniq[static_cast<std::size_t>(k)].val;
		p[uniq[static_cast<std::size_t>(k)].col + 1] += 1;
	}
	for (int c = 0; c < ncol; ++c) {
		p[c + 1] += p[c];
	}

	Rcpp::S4 mat("dgCMatrix");
	mat.slot("i") = i;
	mat.slot("p") = p;
	mat.slot("Dim") = Rcpp::IntegerVector::create(n, ncol);
	mat.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);
	mat.slot("x") = x;
	return mat;
}
