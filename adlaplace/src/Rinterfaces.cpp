#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/runtime/interfaces.hpp"
#include "adlaplace/ompad.hpp"

double jointLogDens(SEXP ad_fun, const Rcpp::NumericVector& x, SEXP Sgroups);
Rcpp::NumericVector grad(SEXP ad_fun, const Rcpp::NumericVector& x, bool inner, SEXP Sgroups);
Rcpp::DataFrame hessian(SEXP ad_fun, const Rcpp::NumericVector& x, bool inner, SEXP Sgroups, const bool verbose);
Rcpp::NumericVector traceHinvT(SEXP ad_fun, const Rcpp::NumericVector& x, const Rcpp::S4& LinvPt, const Rcpp::S4& LinvPtColumns, const int num_threads, SEXP Sgroups);


//' C++ backend entry points
//'
//' Low-level C++ entry points exposed to R via Rcpp.
//' These create and operate on backend state returned by \code{getAdFun_r()}.
//' For the default backend this state is a list with an opaque external pointer
//' plus sparsity/Hessian metadata.
//'
//' @param data An R list containing model data objects required by the backend
//'   (used by \code{getAdFun_r()}).
//' @param config An R list of configuration options required by the backend
//'   (used by \code{getAdFun_r()}).
//' @param ad_fun Backend object returned by \code{getAdFun_r()}.
//'   For the default backend this can be either the full returned list
//'   (containing \code{adFun}, \code{sparsity}, \code{hessians}) or the
//'   external pointer in \code{$adFun}.
//' @param x Numeric parameter vector of length \code{Nparams}.
//' @param inner Logical scalar. If \code{TRUE}, evaluate inner-\eqn{\gamma}
//'   derivatives; otherwise evaluate outer/full derivatives.
//' @param Sgroups Optional integer vector of 0-based group indices to evaluate.
//'   If omitted, uses all groups \code{0:(Ngroups-1)}.
//' @param LinvPt Sparse \code{dgCMatrix} for columns of
//'   \eqn{P^\top L^{-1} D^{-1/2}} (or equivalent) used in trace contractions.
//' @param LinvPtColumns Sparse \code{ngCMatrix}/\code{dgCMatrix} mapping
//'   selected columns of \code{LinvPt} to each group.
//' @param verbose Logical scalar indicating whether to print verbose output.
//' @param num_threads Integer specifying the number of threads to use for parallel computation.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun_r}: backend object. For the default backend: a list
//'     with \code{adFun} (external pointer), \code{sparsity}, and
//'     \code{hessians}.
//'   \item \code{jointLogDens}: scalar objective value summed over groups.
//'   \item \code{grad}: numeric gradient vector.
//'   \item \code{hessian}: data frame with columns group, row, col, value.
//'   \item \code{traceHinvT}: numeric vector of third-derivative contractions.
//' }
//'
//' @details
//' In the default backend, \code{$adFun} is an opaque external pointer and not
//' user-modifiable. It may hold substantial memory (AD tapes and caches).
//' Do not save backend objects across R sessions.
//'
//' @name adlaplace_cpp

//' @rdname adlaplace_cpp
// [[Rcpp::export]]
SEXP getAdFun(
  Rcpp::List data,
  Rcpp::List config)
{
  return getAdFun_h(data, config);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  SEXP Sgroups = R_NilValue) {

  adlaplace_adpack_handle* h = get_handle(ad_fun);

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  double total = 0.0;
  for (size_t g : groups) {
    double fg = 0.0;
    int gi = static_cast<int>(g);
    const int rc = h->api->f(h->ctx, &gi, x.begin(), &fg);
    if (rc != 0) {
      Rcpp::stop("backend api->f failed for group %d with code %d", gi, rc);
    }
    total += fg;
  }
  return total;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  bool inner = false,
  SEXP Sgroups = R_NilValue) {
  adlaplace_adpack_handle* h = get_handle(ad_fun);
  if (!h->api->f_grad) {
    Rcpp::stop("ad_fun api->f_grad is NULL");
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  Rcpp::NumericVector grad_out(Nparams, 0.0);
  double f_total = 0.0;
  const bool inner_local = inner;
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  for (size_t g : groups) {
    int gi = static_cast<int>(g);
    const int rc = h->api->f_grad(
      h->ctx, &gi, x.begin(), &inner_local, &f_total, grad_out.begin()
    );
    if (rc != 0) {
      Rcpp::stop("backend api->f_grad failed for group %d with code %d", gi, rc);
    }
  }
  return grad_out;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame hessian(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  bool inner = false,
  SEXP Sgroups = R_NilValue,
  const bool verbose = false) {
  adlaplace_adpack_handle* h = get_handle(ad_fun);
  if (!h->api->f_grad_hess) {
    Rcpp::stop("ad_fun api->f_grad_hess is NULL");
  }
  if (!h->api->get_sizes) {
    Rcpp::stop("ad_fun api->get_sizes is NULL");
  }
  if (!h->api->get_sparse_sizes) {
    Rcpp::stop("ad_fun api->get_sparse_sizes is NULL");
  }
  if (!h->api->get_sparse_pattern) {
    Rcpp::stop("ad_fun api->get_sparse_pattern is NULL");
  }
  if (verbose) {
    Rcpp::Rcout << "Starting Hessian computation..." << std::endl;
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  const bool inner_local = inner;
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  // First pass: get sparsity info for each group and compute total size
  std::vector<int> group_nnz;
  std::vector<int*> group_pattern_rows;
  std::vector<int*> group_pattern_cols;
  size_t total_nnz = 0;

  for (size_t g : groups) {
    int gi = static_cast<int>(g);
    
    int n_inner = 0, n_outer = 0;
    int nnz_grad_inner = 0, nnz_grad_outer = 0;
    int nnz_hes_inner = 0, nnz_hes_outer = 0;
    const int rc_group_sizes = h->api->get_sparse_sizes(
      h->ctx, &gi, &n_inner, &n_outer,
      &nnz_grad_inner, &nnz_grad_outer,
      &nnz_hes_inner, &nnz_hes_outer
    );
    if (rc_group_sizes != 0) {
      Rcpp::stop("backend api->get_sparse_sizes failed for group %d with code %d", gi, rc_group_sizes);
    }
    
    int nnz_hes = inner ? nnz_hes_inner : nnz_hes_outer;
    group_nnz.push_back(nnz_hes);
    total_nnz += nnz_hes;
    
    int* pattern_grad_inner = new int[nnz_grad_inner];
    int* pattern_grad_outer = new int[nnz_grad_outer];
    int* pattern_hes_inner_row = new int[nnz_hes_inner];
    int* pattern_hes_inner_col = new int[nnz_hes_inner];
    int* pattern_hes_outer_row = new int[nnz_hes_outer];
    int* pattern_hes_outer_col = new int[nnz_hes_outer];
    
    const int rc_pattern = h->api->get_sparse_pattern(
      h->ctx, &gi,
      pattern_grad_inner, pattern_grad_outer,
      pattern_hes_inner_row, pattern_hes_inner_col,
      pattern_hes_outer_row, pattern_hes_outer_col
    );
    
    // Store the pattern for this group
    group_pattern_rows.push_back(inner ? pattern_hes_inner_row : pattern_hes_outer_row);
    group_pattern_cols.push_back(inner ? pattern_hes_inner_col : pattern_hes_outer_col);
    
    // Clean up arrays we don't need
    delete[] pattern_grad_inner;
    delete[] pattern_grad_outer;
    if (inner) {
      delete[] pattern_hes_outer_row;
      delete[] pattern_hes_outer_col;
    } else {
      delete[] pattern_hes_inner_row;
      delete[] pattern_hes_inner_col;
    }
    
    if (rc_pattern != 0) {
      Rcpp::stop("backend api->get_sparse_pattern failed for group %d with code %d", gi, rc_pattern);
    }
  }

  // Allocate buffer for all Hessian values
  Rcpp::NumericVector hes_values(total_nnz, 0.0);
  Rcpp::NumericVector grad_scratch(Nparams, 0.0);
  double f_total = 0.0;

  // Second pass: compute Hessian for each group
  size_t offset = 0;
  std::vector<int> group;
  std::vector<int> row_idx;
  std::vector<int> col_idx;

  for (size_t g_idx = 0; g_idx < groups.size(); g_idx++) {
    size_t g = groups[g_idx];
    int gi = static_cast<int>(g);
    int nnz_hes = group_nnz[g_idx];
    
    // Create map: local index -> position in hes_values
    int* map = new int[nnz_hes];
    for (int i = 0; i < nnz_hes; i++) {
      map[i] = static_cast<int>(offset + i);
    }
    
    const int rc = h->api->f_grad_hess(
      h->ctx,
      &gi,
      x.begin(),
      &inner_local,
      &f_total,
      grad_scratch.begin(),
      hes_values.begin(),
      map
    );
    delete[] map;
    
    if (rc != 0) {
      Rcpp::stop("backend api->f_grad_hess failed for group %d with code %d", gi, rc);
    }
    
    // Store group and indices
    int* rows = group_pattern_rows[g_idx];
    int* cols = group_pattern_cols[g_idx];
    for (int i = 0; i < nnz_hes; i++) {
      group.push_back(static_cast<int>(g));
      row_idx.push_back(rows[i]);
      col_idx.push_back(cols[i]);
    }
    
    // Clean up pattern arrays
    delete[] rows;
    delete[] cols;
    
    offset += nnz_hes;
  }

  if (verbose) {
    Rcpp::Rcout << "Hessian computation completed successfully" << std::endl;
  }

  // Build data frame
  Rcpp::DataFrame result = Rcpp::DataFrame::create(
    Rcpp::Named("group") = group,
    Rcpp::Named("row") = row_idx,
    Rcpp::Named("col") = col_idx,
    Rcpp::Named("value") = hes_values
  );

  return result;

}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT( // to do: pass num threads
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  const int num_threads,
  SEXP Sgroups
) {
  adlaplace_adpack_handle* h = get_handle(ad_fun);
  if (!h->api->trace_hinv_t) {
    Rcpp::stop("ad_fun api->trace_hinv_t is NULL");
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  Rcpp::IntegerVector LinvPt_p = LinvPt.slot("p");
  Rcpp::IntegerVector LinvPt_i = LinvPt.slot("i");
  Rcpp::NumericVector LinvPt_x = LinvPt.slot("x");
  Rcpp::IntegerVector LinvPt_Dim = LinvPt.slot("Dim");

  Rcpp::IntegerVector LinvPtColumns_p = LinvPtColumns.slot("p");
  Rcpp::IntegerVector LinvPtColumns_i = LinvPtColumns.slot("i");

  const size_t LinvPt_ncol = static_cast<size_t>(LinvPt_Dim[1]);
  const size_t LinvPt_p_len = static_cast<size_t>(LinvPt_p.size());
  const size_t LinvPt_i_len = static_cast<size_t>(LinvPt_i.size());
  const size_t LinvPt_x_len = static_cast<size_t>(LinvPt_x.size());
  const size_t LinvPtColumns_p_len = static_cast<size_t>(LinvPtColumns_p.size());
  const size_t LinvPtColumns_i_len = static_cast<size_t>(LinvPtColumns_i.size());

  std::vector<double> trace_accum(Nparams, 0.0);
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  int rc_error = 0;
  int rc_group = -1;
  const int n_groups = static_cast<int>(groups.size());

  if (n_groups > 0) {
    cppad_parallel_setup(static_cast<std::size_t>(num_threads));

#pragma omp parallel num_threads(num_threads)
    {
      std::vector<double> trace_local(Nparams, 0.0);
      int rc_local = 0;
      int rc_group_local = -1;

#pragma omp for schedule(static,1)
      for (int gidx = 0; gidx < n_groups; ++gidx) {
        const int gi = static_cast<int>(groups[static_cast<std::size_t>(gidx)]);
        const int rc = h->api->trace_hinv_t(
          h->ctx,
          &gi,
          x.begin(),
          LinvPt_p.begin(),
          LinvPt_i.begin(),
          LinvPt_x.begin(),
          LinvPt_ncol,
          LinvPt_p_len,
          LinvPt_i_len,
          LinvPt_x_len,
          LinvPtColumns_p.begin(),
          LinvPtColumns_i.begin(),
          LinvPtColumns_p_len,
          LinvPtColumns_i_len,
          trace_local.data()
        );
        if (rc != 0 && rc_local == 0) {
          rc_local = rc;
          rc_group_local = gi;
        }
      }

#pragma omp critical
      {
        if (rc_error == 0 && rc_local != 0) {
          rc_error = rc_local;
          rc_group = rc_group_local;
        }
        for (size_t d = 0; d < Nparams; ++d) {
          trace_accum[d] += trace_local[d];
        }
      }
    }
  }

  if (rc_error != 0) {
    Rcpp::stop("backend api->trace_hinv_t failed for group %d with code %d", rc_group, rc_error);
  }

  Rcpp::NumericVector trace_out(Nparams);
  for (size_t d = 0; d < Nparams; ++d) {
    trace_out[static_cast<R_xlen_t>(d)] = trace_accum[d];
  }

  return trace_out;
}
