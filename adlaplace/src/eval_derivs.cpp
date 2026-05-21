#include <Rcpp.h>

#include <map>
#include <vector>

#include "adlaplace/runtime/interfaces.hpp"

//' @title C++ backend entry points
//' @name adlaplace_cpp
//' @description Low-level C++ entry points exposed to R via Rcpp.
//'
//' @section Sign convention:
//' \code{jointLogDens()}, \code{grad()}, and \code{hessian()} evaluate the
//' \strong{joint log density} \eqn{\ell(x)} and its derivatives (maximization sign).
//' \code{all_derivs()} and the objective inside \code{inner_opt()} use the
//' \strong{negative} log density \eqn{-\ell(x)} and its derivatives for
//' \pkg{trustOptim} minimization. At the same \code{x} and \code{ad_fun},
//' \code{all_derivs()$fval == -jointLogDens(ad_fun, x)},
//' \code{all_derivs()$gradient == -grad(ad_fun, x)}, and
//' \code{all_derivs()$hessian == -hessian(ad_fun, x)} (outer, full parameter vector).
//' @param ad_fun External pointer or list from \code{getAdFun()}.
//' @param x Numeric parameter vector of length \code{Nparams}.
//' @param Sgroups Optional integer vector of 0-based group indices.
//' @param inner Logical scalar for inner-\eqn{\gamma} vs outer derivatives.
//' @param verbose Logical passed to \code{hessian()}.
//' @param LinvPt,LinvPtColumns,num_threads See \code{traceHinvT()}.
//'
//' @rdname adlaplace_cpp
//' @return Scalar **joint log density** \eqn{\ell(x)} (sum over shards; not negative log density).
//' @export
// [[Rcpp::export]]
double jointLogDens(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  SEXP Sgroups = R_NilValue) {

  ad_groups* groups = resolve_ad_groups(ad_fun);
  const size_t Ngroups = groups->fun.size();
  if (Ngroups == 0) Rcpp::stop("ad_groups.fun is empty");
  const size_t Nparams = pack_ctx(groups->fun[0]->ctx)->x.size();
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> group_idx = resolve_groups(Ngroups, Sgroups_vec);

  double total = 0.0;
  for (size_t g : group_idx) {
    adlaplace_adpack_handle* h = shard_handle(groups, g);
    double fg = 0.0;
    const int rc = h->api->f(h->ctx, x.begin(), &fg);
    if (rc != 0) {
      Rcpp::stop("backend api->f failed for group %d with code %d", (int)g, rc);
    }
    total += fg;
  }
  return total;
}

//' @rdname adlaplace_cpp
//' @return Gradient of **log density** \eqn{\nabla \ell(x)} (not negative log density).
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  SEXP Sgroups = R_NilValue,
  bool inner = false) {

  ad_groups* groups = resolve_ad_groups(ad_fun);
  const size_t Nparams = x.size();
  const size_t Ngroups = groups->fun.size();
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> group_idx = resolve_groups(Ngroups, Sgroups_vec);

  Rcpp::NumericVector grad_out(Nparams, 0.0);
  double f_dummy = 0.0;

  for (size_t g : group_idx) {
    adlaplace_adpack_handle* h = shard_handle(groups, g);
    if (h->api->f_grad(h->ctx, x.begin(), &inner, &f_dummy, grad_out.begin()) != 0)
      Rcpp::stop("backend api->f_grad failed for group %d", (int)g);
  }
  return grad_out;
}

//' @rdname adlaplace_cpp
//' @return Sparse Hessian of **log density** \eqn{\nabla^2 \ell(x)} (not negative log density).
//' @export
// [[Rcpp::export]]
Rcpp::S4 hessian(
  SEXP ad_fun,
  const Rcpp::NumericVector& x,
  SEXP Sgroups = R_NilValue,
  bool inner = false,
  const bool verbose = false) {
  ad_groups* groups = resolve_ad_groups(ad_fun);
  if (groups->fun.empty() || !groups->fun[0]->api->f_grad_hess) {
    Rcpp::stop("ad_fun api->f_grad_hess is NULL");
  }

  const size_t Nparams = x.size();
  const size_t Ngroups = groups->fun.size();
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> group_idx = resolve_groups(Ngroups, Sgroups_vec);

  std::map<std::pair<int, int>, double> acc_map;

  if (verbose) Rcpp::Rcout << "Starting Hessian computation..." << std::endl;

  for (size_t g : group_idx) {
    adlaplace_adpack_handle* h = shard_handle(groups, g);
    int n_inner, n_outer, nnz_grad_i, nnz_grad_o, nnz_hes_i, nnz_hes_o;
    if (h->api->get_sparse_sizes(h->ctx, &n_inner, &n_outer, &nnz_grad_i, &nnz_grad_o, &nnz_hes_i, &nnz_hes_o) != 0)
      Rcpp::stop("backend api->get_sparse_sizes failed for group %d", (int)g);

    int nnz = inner ? nnz_hes_i : nnz_hes_o;
    if (nnz == 0) continue;

    std::vector<int> p_gi(nnz_grad_i), p_go(nnz_grad_o);
    std::vector<int> p_hir(nnz_hes_i), p_hic(nnz_hes_i);
    std::vector<int> p_hor(nnz_hes_o), p_hoc(nnz_hes_o);

    if (h->api->get_sparse_pattern(h->ctx,
        p_gi.data(), p_go.data(),
        p_hir.data(), p_hic.data(),
        p_hor.data(), p_hoc.data()) != 0)
      Rcpp::stop("backend api->get_sparse_pattern failed for group %d", (int)g);

    const int* rows = inner ? p_hir.data() : p_hor.data();
    const int* cols = inner ? p_hic.data() : p_hoc.data();

    std::vector<double> vals(nnz);
    std::vector<int> map(nnz);
    for (int i = 0; i < nnz; ++i) map[i] = i;

    std::vector<double> grad_scratch(Nparams, 0.0);
    double f_dummy = 0.0;
    if (h->api->f_grad_hess(h->ctx, x.begin(), &inner, &f_dummy, grad_scratch.data(), vals.data(), map.data()) != 0)
      Rcpp::stop("backend api->f_grad_hess failed for group %d", (int)g);

    for (int i = 0; i < nnz; ++i) {
      acc_map[{rows[i], cols[i]}] += vals[i];
    }
  }

  if (verbose) Rcpp::Rcout << "Hessian computation completed successfully" << std::endl;

  std::vector<int> agg_row, agg_col;
  std::vector<double> agg_value;
  agg_row.reserve(acc_map.size());
  agg_col.reserve(acc_map.size());
  agg_value.reserve(acc_map.size());

  for (const auto& kv : acc_map) {
    agg_row.push_back(kv.first.first);
    agg_col.push_back(kv.first.second);
    agg_value.push_back(kv.second);
  }

  static Rcpp::Function sparseMatrix =
    Rcpp::Environment::namespace_env("Matrix")["sparseMatrix"];
  const int n = static_cast<int>(Nparams);

  return sparseMatrix(
    Rcpp::Named("i") = Rcpp::wrap(agg_row),
    Rcpp::Named("j") = Rcpp::wrap(agg_col),
    Rcpp::Named("x") = Rcpp::wrap(agg_value),
    Rcpp::Named("index1") = false,
    Rcpp::Named("symmetric") = true,
    Rcpp::Named("dims") = Rcpp::IntegerVector::create(n, n)
  );
}
