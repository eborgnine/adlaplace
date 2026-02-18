#include <Rcpp.h>

#include <cstddef>
#include <algorithm>
#include <limits>
#include <vector>

namespace {

inline void require_name(const Rcpp::List& x, const char* name) {
  if (!x.containsElementNamed(name)) {
    Rcpp::stop("sparsity shard missing `%s`", name);
  }
}

struct SparseTriplet {
  int shard;
  int row;
  int col;
  int local;
  long long cell;
};

inline std::vector<long long> sorted_unique_cells(const std::vector<SparseTriplet>& x) {
  std::vector<long long> cells;
  cells.reserve(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) cells.push_back(x[k].cell);
  std::sort(cells.begin(), cells.end());
  cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
  return cells;
}

inline int cell_sparse_id(const std::vector<long long>& sorted_cells, long long cell) {
  const auto it = std::lower_bound(sorted_cells.begin(), sorted_cells.end(), cell);
  if (it == sorted_cells.end() || *it != cell) {
    Rcpp::stop("internal error: cell not found");
  }
  return static_cast<int>(it - sorted_cells.begin());
}

}  // namespace

// cpp version of hessianMap, written by codex.

// [[Rcpp::export]]
Rcpp::List hessianMapC(
  const Rcpp::List& sparsity_list,
  const int Nbeta,
  const int Ngamma,
  const int Ntheta
) {
  if (Nbeta < 0 || Ngamma < 0 || Ntheta < 0) {
    Rcpp::stop("Nbeta, Ngamma, and Ntheta must be non-negative");
  }

  const int Ngroups = static_cast<int>(sparsity_list.size());
  const int Nparams = Nbeta + Ngamma + Ntheta;

  std::vector<SparseTriplet> outer_entries;
  std::vector<SparseTriplet> inner_entries;
  outer_entries.reserve(1024);
  inner_entries.reserve(1024);

  for (int shard = 0; shard < Ngroups; ++shard) {
    const Rcpp::List shard_list = Rcpp::as<Rcpp::List>(sparsity_list[shard]);

    require_name(shard_list, "row_hess");
    require_name(shard_list, "col_hess");
    require_name(shard_list, "row_hess_inner");
    require_name(shard_list, "col_hess_inner");

    const Rcpp::IntegerVector row_hess = shard_list["row_hess"];
    const Rcpp::IntegerVector col_hess = shard_list["col_hess"];
    const R_xlen_t nouter = col_hess.size();
    if (row_hess.size() != nouter) {
      Rcpp::stop("row_hess/col_hess length mismatch for shard %d", shard);
    }

    for (R_xlen_t k = 0; k < nouter; ++k) {
      const int row = row_hess[k];
      const int col = col_hess[k];
      const long long cell =
        static_cast<long long>(row) + static_cast<long long>(col) * static_cast<long long>(Nparams);

      outer_entries.push_back(SparseTriplet{
        shard, row, col, static_cast<int>(k), cell
      });
    }

    const Rcpp::IntegerVector row_hess_inner = shard_list["row_hess_inner"];
    const Rcpp::IntegerVector col_hess_inner = shard_list["col_hess_inner"];
    const R_xlen_t ninner = col_hess_inner.size();
    if (row_hess_inner.size() != ninner) {
      Rcpp::stop("row_hess_inner/col_hess_inner length mismatch for shard %d", shard);
    }

    for (R_xlen_t k = 0; k < ninner; ++k) {
      const int row = row_hess_inner[k];
      const int col = col_hess_inner[k];
      const int row_inner = row - Nbeta;
      const int col_inner = col - Nbeta;

      if (row_inner < 0 || row_inner >= Ngamma || col_inner < 0 || col_inner >= Ngamma) {
        Rcpp::stop("inner hessian index out of gamma range in shard %d", shard);
      }

      const long long cell = static_cast<long long>(row_inner) +
        static_cast<long long>(Ngamma) * static_cast<long long>(col_inner);

      inner_entries.push_back(SparseTriplet{
        shard, row_inner, col_inner, static_cast<int>(k), cell
      });
    }
  }

  const std::vector<long long> outer_cells = sorted_unique_cells(outer_entries);
  const std::vector<long long> inner_cells = sorted_unique_cells(inner_entries);

  std::vector<int> outer_shard;
  std::vector<int> outer_cell_sparse;
  std::vector<int> outer_local;
  outer_shard.reserve(outer_entries.size());
  outer_cell_sparse.reserve(outer_entries.size());
  outer_local.reserve(outer_entries.size());

  std::vector<int> inner_shard;
  std::vector<int> inner_cell_sparse;
  std::vector<int> inner_local;
  inner_shard.reserve(inner_entries.size());
  inner_cell_sparse.reserve(inner_entries.size());
  inner_local.reserve(inner_entries.size());

  for (std::size_t k = 0; k < outer_entries.size(); ++k) {
    outer_shard.push_back(outer_entries[k].shard);
    outer_cell_sparse.push_back(cell_sparse_id(outer_cells, outer_entries[k].cell));
    outer_local.push_back(outer_entries[k].local);
  }
  for (std::size_t k = 0; k < inner_entries.size(); ++k) {
    inner_shard.push_back(inner_entries[k].shard);
    inner_cell_sparse.push_back(cell_sparse_id(inner_cells, inner_entries[k].cell));
    inner_local.push_back(inner_entries[k].local);
  }

  std::vector<int> outer_unique_row;
  std::vector<int> outer_unique_col;
  std::vector<int> outer_unique_x;
  std::vector<int> outer_seen(outer_cells.size(), 0);
  for (std::size_t k = 0; k < outer_entries.size(); ++k) {
    const int sid = outer_cell_sparse[k];
    if (!outer_seen[sid]) {
      outer_seen[sid] = 1;
      outer_unique_row.push_back(outer_entries[k].row);
      outer_unique_col.push_back(outer_entries[k].col);
      outer_unique_x.push_back(sid);
    }
  }

  std::vector<int> inner_unique_row;
  std::vector<int> inner_unique_col;
  std::vector<int> inner_unique_x;
  std::vector<int> inner_seen(inner_cells.size(), 0);
  for (std::size_t k = 0; k < inner_entries.size(); ++k) {
    const int sid = inner_cell_sparse[k];
    if (!inner_seen[sid]) {
      inner_seen[sid] = 1;
      inner_unique_row.push_back(inner_entries[k].row);
      inner_unique_col.push_back(inner_entries[k].col);
      inner_unique_x.push_back(sid);
    }
  }

  Rcpp::Function sparse_matrix = Rcpp::Environment::namespace_env("Matrix")["sparseMatrix"];
  const bool index1 = false;
  const bool symmetric = true;

  const Rcpp::S4 hessian = Rcpp::as<Rcpp::S4>(sparse_matrix(
    Rcpp::Named("i") = Rcpp::wrap(outer_unique_row),
    Rcpp::Named("j") = Rcpp::wrap(outer_unique_col),
    Rcpp::Named("x") = Rcpp::wrap(outer_unique_x),
    Rcpp::Named("symmetric") = symmetric,
    Rcpp::Named("index1") = index1,
    Rcpp::Named("dims") = Rcpp::IntegerVector::create(Nparams, Nparams)
  ));

  const Rcpp::S4 hessian_inner = Rcpp::as<Rcpp::S4>(sparse_matrix(
    Rcpp::Named("i") = Rcpp::wrap(inner_unique_row),
    Rcpp::Named("j") = Rcpp::wrap(inner_unique_col),
    Rcpp::Named("x") = Rcpp::wrap(inner_unique_x),
    Rcpp::Named("symmetric") = symmetric,
    Rcpp::Named("index1") = index1,
    Rcpp::Named("dims") = Rcpp::IntegerVector::create(Ngamma, Ngamma)
  ));

  const Rcpp::S4 hessian_outer_map = Rcpp::as<Rcpp::S4>(sparse_matrix(
    Rcpp::Named("i") = Rcpp::wrap(outer_cell_sparse),
    Rcpp::Named("j") = Rcpp::wrap(outer_shard),
    Rcpp::Named("x") = Rcpp::wrap(outer_local),
    Rcpp::Named("index1") = index1,
    Rcpp::Named("dims") = Rcpp::IntegerVector::create(
      static_cast<int>(outer_cells.size()), Ngroups
    )
  ));

  const Rcpp::S4 hessian_inner_map = Rcpp::as<Rcpp::S4>(sparse_matrix(
    Rcpp::Named("i") = Rcpp::wrap(inner_cell_sparse),
    Rcpp::Named("j") = Rcpp::wrap(inner_shard),
    Rcpp::Named("x") = Rcpp::wrap(inner_local),
    Rcpp::Named("index1") = index1,
    Rcpp::Named("dims") = Rcpp::IntegerVector::create(
      static_cast<int>(inner_cells.size()), Ngroups
    )
  ));

  const Rcpp::IntegerVector inner_map_i = hessian_inner_map.slot("i");
  const Rcpp::NumericVector inner_hessian_x_num = hessian_inner.slot("x");
  Rcpp::IntegerVector inner_hessian_x(inner_hessian_x_num.size());
  for (R_xlen_t i = 0; i < inner_hessian_x_num.size(); ++i) {
    inner_hessian_x[i] = static_cast<int>(inner_hessian_x_num[i]);
  }
  Rcpp::Function match_fn("match");
  Rcpp::IntegerVector inner_global = match_fn(inner_map_i, inner_hessian_x);
  for (R_xlen_t i = 0; i < inner_global.size(); ++i) {
    inner_global[i] = inner_global[i] - 1;
  }

  const Rcpp::List result_map = Rcpp::List::create(
    Rcpp::Named("outer") = Rcpp::List::create(
      Rcpp::Named("p") = Rcpp::as<Rcpp::IntegerVector>(hessian_outer_map.slot("p")),
      Rcpp::Named("local") = Rcpp::as<Rcpp::IntegerVector>(hessian_outer_map.slot("x")),
      Rcpp::Named("global") = Rcpp::as<Rcpp::IntegerVector>(hessian_outer_map.slot("i")),
      Rcpp::Named("dims") = Rcpp::as<Rcpp::IntegerVector>(hessian_outer_map.slot("Dim"))
    ),
    Rcpp::Named("inner") = Rcpp::List::create(
      Rcpp::Named("p") = Rcpp::as<Rcpp::IntegerVector>(hessian_inner_map.slot("p")),
      Rcpp::Named("local") = Rcpp::as<Rcpp::IntegerVector>(hessian_inner_map.slot("x")),
      Rcpp::Named("global") = inner_global,
      Rcpp::Named("dims") = Rcpp::IntegerVector::create(
        Rf_length(hessian_inner.slot("x")),
        Rcpp::as<Rcpp::IntegerVector>(hessian_inner_map.slot("Dim"))[1]
      )
    )
  );

  return Rcpp::List::create(
    Rcpp::Named("hessian") = Rcpp::List::create(
      Rcpp::Named("outer") = hessian,
      Rcpp::Named("inner") = hessian_inner
    ),
    Rcpp::Named("map") = result_map
  );
}
