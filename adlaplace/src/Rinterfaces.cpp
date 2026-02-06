#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/creators/R_interfaces.hpp"






//' C++ backend entry points
//'
//' Low-level C++ entry points exposed to R via Rcpp. These create and operate on
//' an opaque AD “pack” (external pointer) used to evaluate the objective,
//' sparse gradient values, and sparse Hessian values for a selected group.
//'
//' Indices passed in `pattern`, `row`, and `col` are **0-based** and must be
//' sorted with no duplicates (for gradients) or paired consistently (for Hessians).
//'
//' @param data An R list containing model data objects required by the backend.
//' @param config An R list of configuration options required by the backend.
//' @param parameters Numeric vector of parameters.
//' @param i Integer index selecting which group/tape to evaluate.
//' @param adPack External pointer (`externalptr`) returned by \code{getAdFun()}.
//' @param pattern Integer vector of 0-based column indices specifying which
//'   gradient entries to compute (sorted, unique).
//' @param row Integer vector of 0-based row indices for Hessian entries.
//' @param col Integer vector of 0-based column indices for Hessian entries;
//'   must have the same length as \code{row}.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun}: a list containing an opaque external pointer and
//'     associated metadata.
//'   \item \code{jointLogDens}: a scalar objective value for group \code{i}.
//'   \item \code{grad}: numeric vector of length \code{length(pattern)} with
//'     gradient values in the same order as \code{pattern}.
//'   \item \code{hess}: numeric vector of length \code{length(row)} with Hessian
//'     values in the same order as the \code{(row, col)} pairs.
//' }
//'
//' @details
//' The external pointer returned by \code{getAdFun} is opaque and not
//' user-modifiable. It may hold substantial memory (AD tapes, sparsity patterns,
//' work caches). Do not save it across R sessions.
//'
//' @name adlaplace_cpp

// in external_create.hpp
SEXP getAdFun_h(
	const Data& data,
	const Config& config);

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
SEXP getAdFun(
	Rcpp::List data, 
	Rcpp::List config)
{

	SEXP result = getAdFun_h(data, config);

	return result;
}



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(const Rcpp::NumericVector& x, SEXP backendContext) {

  adlaplace_adpack_handle* h = get_handle(backendContext);

  const size_t Ngroups = h->Ngroups;
  const size_t Nparams = h->Nparams;

  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  // Evaluate each group's contribution
  for (size_t g = 0; g < Ngroups; ++g) {
    double fg = 0.0;
    int gi = static_cast<int>(g);
    int rc = h->api->f(h->ctx, &gi, x.begin(), &fg);
    if (rc != 0) {
      Rcpp::stop("backend api->f failed for group %d with code %d", gi, rc);
    }
  }

  return fg;
}
