#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/adpack_api.h"
#include "adlaplace/adpack_handle.h"

#include "adlaplace/utils.hpp"


inline CPPAD_TESTVECTOR(double) as_cppad_vector(
  const Rcpp::NumericVector& v
  ) {
  const size_t n = v.size();
  CPPAD_TESTVECTOR(double) out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    out[i] = v[i];
  }
  return out;
}


#ifdef UNDEF

void grad(
	const CPPAD_TESTVECTOR(double)& parameters, 
	const Config& config,
	SEXP adPack,
	CPPAD_TESTVECTOR(double) &result,
	const bool inner=false) {


	const size_t Nout = inner?config.Ngamma:config.Nparams;
	if(result.size() != Nout) {
		Rcpp::Rcout << "result wrong size " << result.size() << " inner " << inner << " Ngamma " << config.Ngamma <<
		" Nparams " << config.Nparams << "\n";
	}

	CPPAD_TESTVECTOR(double) resultLocal(result.size());

	for(size_t Dgroup: config.Sgroups) {

		const CscPattern& pattern = inner?config.match.grad_inner:config.match.grad;
		const size_t Nhere = pattern.p[Dgroup+1] - pattern.p[Dgroup];
		resultLocal.resize(Nhere);
		grad_api(parameters, Dgroup, adPack, inner, resultLocal);
		for(size_t Dlocal=0,Di=pattern.p[Dgroup];Dlocal < Nhere; Dlocal++,Di++) {

			result[pattern.i[Di]] += resultLocal[Dlocal];

#ifdef DEBUG
			if (pattern.i[Di] >= config.Nparams) {
				Rcpp::stop("grad: match index out of range");
			};
#endif

		}
	}
}

void hess(
	const CPPAD_TESTVECTOR(double)& parameters,
	const Config& config,
	SEXP adPack,
	CPPAD_TESTVECTOR(double) &result,
	const bool inner = false
	) {

	const size_t Nout = inner?config.hessian_inner.nnz():config.hessian.nnz();
	if(result.size() != Nout) {
		Rcpp::Rcout << "reuslt wrong size " << result.size() << " inner " << inner << " innernnz " << config.hessian_inner.nnz() <<
		" hessiannnz " << config.hessian.nnz() << "\n";
	}

	CPPAD_TESTVECTOR(double) resultLocal(result.size());
	CPPAD_TESTVECTOR(double) result_grad(0);


	const CscPattern& pattern = inner?config.match.hessian_inner:config.match.hessian;

	for(size_t Dgroup: config.Sgroups) {

		const std::size_t start = pattern.p[Dgroup];
		const std::size_t end   = pattern.p[Dgroup + 1];
		const std::size_t Nhere = end - start;

		resultLocal.resize(Nhere);

    // Fill group-local Hessian values
		hess_api(parameters, Dgroup, adPack, inner, result_grad, resultLocal);

		for(size_t Di=start;Di < end; Di++) {
			result[pattern.i[Di]] += resultLocal[pattern.x[Di]];
		}
	}
}

#endif


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
Rcpp::List getAdFun_h(
	const Data& data,
	const Config& config);

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List getAdFun(
	Rcpp::List data, 
	Rcpp::List config)
{

	Data dataC(data);
	Config configC(config);

	Rcpp::List result = getAdFun_h(dataC, configC);

	return result;
}

// helper: get handle safely
static inline adlaplace_adpack_handle* get_handle(SEXP handle_sexp) {
  auto* h = static_cast<adlaplace_adpack_handle*>(R_ExternalPtrAddr(handle_sexp));
  if (!h) Rcpp::stop("backendContext handle is NULL (external pointer cleared?)");
  if (!h->api) Rcpp::stop("backendContext handle has NULL api");
  if (!h->ctx) Rcpp::stop("backendContext handle has NULL ctx");
  if (!h->api->f) Rcpp::stop("backendContext api->f is NULL");
  return h;
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
