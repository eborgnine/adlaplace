#include<Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/adfun.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List sparsity_grad(
	SEXP adPack,
	const Rcpp::NumericVector parameters
	) {
	
	const size_t Nparams = parameters.size();
	CppAD::vector<double> parametersC(Nparams);
	for (size_t j = 0; j < Nparams; ++j)
		parametersC[j] = parameters[j];

	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	std::vector<GroupPack>& adPackC = *xp;

	const size_t Nfun = adPackC.size();
	Rcpp::List result(Nfun);

	CppAD::sparse_rc< CppAD::vector<size_t> > pattern_in;
    pattern_in.resize(Nparams, Nparams, Nparams);   // n nonzeros
    for (size_t j = 0; j < Nparams; ++j) {
      pattern_in.set(j, j, j);   // (row=j, col=j, index=j)
  }

    bool transpose    = false;   // want pattern_out as m × n
    bool dependency   = false;   // standard Jacobian sparsity (not dependency)
    bool internal_bool = false;  // let CppAD choose representation; updated on return

	CPPAD_TESTVECTOR(bool) select_domain(Nparams), select_range(1);
		for (size_t j = 0; j < Nparams; ++j) {
			select_domain[j] = true;
		}
    select_range[0] = true; // scalar output


    for(size_t D=0;D<Nfun;D++) {
    	CppAD::sparse_rc< CppAD::vector<size_t> > pattern_out_grad;
	    CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t) > pattern_out_hess;

    	adPackC[D].fun.for_jac_sparsity(
    		pattern_in,
    		transpose,
    		dependency,
    		internal_bool,
    		pattern_out_grad);

    	adPackC[D].fun.for_hes_sparsity(
    	select_domain,
    	select_range,
    	internal_bool,
    	pattern_out_hess
    	);


    	const CppAD::vector<size_t>& col_grad = pattern_out_grad.col();
    	const size_t Kgrad = pattern_out_grad.nnz();
    	Rcpp::IntegerVector resultHereGrad(Kgrad);

    const CppAD::vector<size_t>& row = pattern_out_hess.row();
    const CppAD::vector<size_t>& col = pattern_out_hess.col();
    const size_t K = pattern_out_hess.nnz();
    	Rcpp::IntegerVector resultHereRow(K);
    	Rcpp::IntegerVector resultHereCol(K);

    	for (size_t k = 0; k < Kgrad; ++k) {
    		resultHereGrad[k] = col_grad[k];
    	}
    	for (size_t k = 0; k < K; ++k) {
    		resultHereRow[k] = row[k];
    		resultHereCol[k] = col[k];
    	}

    	result[D] = Rcpp::List::create(
    		Rcpp::Named("grad") = resultHereGrad,
    		Rcpp::Named("row") = resultHereRow,
    		Rcpp::Named("col") = resultHereCol
    		);
    }
    return(result);
}


//' @export
// [[Rcpp::export]]
SEXP addSparsityToAdFun(
	SEXP adPack,
	Rcpp::List sparsity)
{
	Rcpp::XPtr<std::vector<GroupPack>> xp(adPack);
	std::vector<GroupPack>& adPackC = *xp;

	const size_t Nfun = adPackC.size();
	for(size_t D=0;D<Nfun;D++) {
		
	}

}