#include "adlaplace/sparsity.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List sparsity(
	SEXP adPack,
	const Rcpp::NumericVector parameters,
	const bool verbose = false
	) {
	
	const size_t Nparams = parameters.size();
	CppAD::vector<double> parametersC(Nparams);
	for (size_t j = 0; j < Nparams; ++j)
		parametersC[j] = parameters[j];

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);

	Rcpp::List result = sparsity(*xp, parametersC, verbose);

	return(result);
}

