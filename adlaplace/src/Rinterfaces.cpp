#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "adlaplace/data.hpp"


Rcpp::List getAdFun_api(
	const Data& data,
	const Config& config);

void fval_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const int i, 
	SEXP adPack, 
	double& result);

void grad_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const int i, 
	SEXP adPack,
	bool inner,
	CPPAD_TESTVECTOR(double) &result 
	);

void hess_api(
	const CPPAD_TESTVECTOR(double) &x, 
	const int i, 
	SEXP adPack, 
	bool inner,
	CPPAD_TESTVECTOR(double) &result );

CPPAD_TESTVECTOR(double) thirdDirection_api(
	const CPPAD_TESTVECTOR(double)&  x,
	const CPPAD_TESTVECTOR(double)&  direction,
	const CPPAD_TESTVECTOR(double)&  direction2,  // all zeros
	const CPPAD_TESTVECTOR(double)&  w, // {0.0, 0.0, 1.0}
	const int i, 
	SEXP adPack
	);

inline CPPAD_TESTVECTOR(double)
as_cppad_vector(const Rcpp::NumericVector& v)
{
	const size_t n = v.size();
	CPPAD_TESTVECTOR(double) out(n);
	for (R_xlen_t i = 0; i < n; ++i) {
		out[i] = v[i];
	}
	return out;
}



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



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List getAdFun(
	Rcpp::List data, 
	Rcpp::List config)
{

	Data dataC(data);
	Config configC(config);

	Rcpp::List result = getAdFun_api(dataC, configC);

	return result;
}

double jointLogDens(
	const CPPAD_TESTVECTOR(double) parameters, 
	const Config config,
	SEXP adPack) {

	double result=0.0;
	for(size_t D=0; D<config.Ngroups;D++) {
		double resultHere;
		fval_api(parameters, D, adPack, resultHere);
		result += resultHere;
	}
	return result;
}

CPPAD_TESTVECTOR(double) grad(
	const CPPAD_TESTVECTOR(double)& parameters, 
	const Config& config,
	SEXP adPack,
	const bool inner=false) {

	const size_t Nout = inner?config.Ngamma:config.Nparams;
	CPPAD_TESTVECTOR(double) result(Nout, 0.0);
	CPPAD_TESTVECTOR(double) resultLocal(result.size());

	for(size_t Dgroup=0; Dgroup<config.Ngroups;Dgroup++) {
		const CscPattern& pattern = inner?config.match.grad_inner:config.match.grad;
		const size_t Nhere = pattern.p[Dgroup+1] - pattern.p[Dgroup];
		resultLocal.resize(Nhere);
		grad_api(parameters, Dgroup, adPack, inner, resultLocal);
		for(size_t Dlocal=0,Di=pattern.p[Dgroup];Dlocal < Nhere; Dlocal++,Di++) {
#ifdef DEBUG
			if (pattern.i[Di] >= config.Nparams) {
				Rcpp::stop("grad: match index out of range");
			};
#endif
			result[pattern.i[Di]] += resultLocal[Dlocal];
		}
	}
	return result;
}

CPPAD_TESTVECTOR(double) hess(
	const CPPAD_TESTVECTOR(double)& parameters,
	const Config& config,
	SEXP adPack,
	const bool inner = false
	) {

	const size_t Nout = inner?config.hessian_inner.nnz():config.hessian.nnz();
	CPPAD_TESTVECTOR(double) result(Nout, 0.0);
	CPPAD_TESTVECTOR(double) resultLocal(result.size());


	for (std::size_t Dgroup = 0; Dgroup < config.Ngroups; ++Dgroup) {
		const CscPattern& pattern = inner?config.match.hessian_inner:config.match.hessian;

		const std::size_t start = pattern.p[Dgroup];
		const std::size_t end   = pattern.p[Dgroup + 1];
		const std::size_t Nhere = end - start;

		resultLocal.resize(Nhere);

    // Fill group-local Hessian values
		hess_api(parameters, Dgroup, adPack, inner, resultLocal);

		for(size_t Dlocal=0,Di=start;Dlocal < Nhere; Dlocal++,Di++) {
			result[pattern.i[Di]] += resultLocal[Dlocal];
		}
	}

	return result;
}
//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(
	const Rcpp::NumericVector parameters, 
	const Rcpp::List config,
	SEXP adPack)
{

	const auto parametersC = as_cppad_vector(parameters);
	const Config configC(config);

	double result = jointLogDens(parametersC, configC, adPack);
	return result;
}


//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
	const Rcpp::NumericVector& parameters, 
	const Rcpp::List config,
	SEXP adPack,
	const bool inner)
{

	auto parametersC = as_cppad_vector(parameters);
	Config configC(config);
	const auto resultC = grad(parametersC, configC, adPack, inner);

	return(Rcpp::wrap(resultC));
}



//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::S4 hess(
	const Rcpp::NumericVector& parameters, 
	const Rcpp::List config,
	SEXP adPack,
	const bool inner)
{

	auto parametersC = as_cppad_vector(parameters);
	Config configC(config);

	const auto resultC = hess(parametersC, configC, adPack, inner);

	Rcpp::S4 out = inner ? Rcpp::clone(Rcpp::as<Rcpp::S4>(config["hessian_inner"]))
	: Rcpp::clone(Rcpp::as<Rcpp::S4>(config["hessian"]));


	if (XLENGTH(out.slot("x")) != resultC.size()) {
		Rcpp::stop("hess: template@x length != result length");
	}

	out.slot("x") = Rcpp::wrap(resultC);

	return(out);

}

