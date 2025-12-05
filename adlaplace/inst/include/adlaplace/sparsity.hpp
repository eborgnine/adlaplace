#ifndef SPARSITY_HPP
#define SPARSITY_HPP


//#define DEBUG




inline CppAD::sparse_rc< CppAD::vector<size_t> > build_gradient_pattern_from_R(
  const Rcpp::IntegerVector& grad_index, const size_t N)
{
  const size_t K = static_cast<size_t>(grad_index.size());
  std::set<size_t> idx;

  for (size_t k = 0; k < K; ++k) {
    size_t c = static_cast<size_t>(grad_index[k]);
    idx.insert(c);
  }


  CppAD::sparse_rc< CppAD::vector<size_t> > pattern;
  pattern.resize(1, N, idx.size());
  size_t t = 0;
  for (size_t j : idx) pattern.set(t++, 0, j);

    return pattern;
}


inline CppAD::sparse_rc< CppAD::vector<size_t> > build_hessian_pattern_from_pairs(
  const Rcpp::IntegerVector& row,
  const Rcpp::IntegerVector& col,
  size_t n) {
  const size_t K = row.size();
  CppAD::sparse_rc< CppAD::vector<size_t> > pat;
  pat.resize(n, n, K);
  bool warnRC=true;
  for (size_t k = 0; k < K; ++k) {
    const int r = row[k], c=col[k];
    if(r > c && warnRC) {
      warnRC = false;
      Rcpp::Rcout << "entry lower, should be uppper triangle\n";
    }
    pat.set(k, r, c);
  }
  return pat;
}


inline void addSparsityToAdFun(
  std::vector<GroupPack>& adPack,
  const Rcpp::List& sparsity,
  const bool verbose=false)
{

  const size_t Nfun = adPack.size();
  const size_t Nparams = adPack[0].fun.Domain();

  if(verbose) {
    Rcpp::Rcout << "sparsity " << sparsity.size() << " adpack " << adPack.size() << "\n";
  }

  if(sparsity.size() != Nfun  ) {
    Rcpp::Rcout << "sparsity size smaller than adPack " << sparsity.size() << " " << Nfun << "\n";
  }
  for(size_t D=0;D<Nfun;D++) {

    const Rcpp::List sparseHere = sparsity[D];
    const Rcpp::IntegerVector sparseGrad = sparseHere["grad"];
    const Rcpp::IntegerVector rowOut = sparseHere["i"];
    const Rcpp::IntegerVector colOut = sparseHere["j"];
    const Rcpp::IntegerVector matchHere = sparseHere["match"];
    const Rcpp::IntegerVector pHere = sparseHere["p"];

#ifdef DEBUG
if(verbose) {
    const int sparseGradSize = sparseGrad.size();
    Rcpp::Rcout << "sparsity group " << D;
    Rcpp::Rcout << " grad size " << sparseGradSize << "\n";
}
#endif
    adPack[D].work_grad = CppAD::sparse_jac_work();
    adPack[D].pattern_grad = build_gradient_pattern_from_R(sparseGrad, Nparams);
    adPack[D].out_grad = 
      CppAD::sparse_rcv< CppAD::vector<size_t>, CppAD::vector<double> >(
        adPack[D].pattern_grad); 

    adPack[D].work_hess= CppAD::sparse_hes_work();
        adPack[D].pattern_hess=build_hessian_pattern_from_pairs(rowOut, colOut, Nparams);
        adPack[D].out_hess = 
          CppAD::sparse_rcv< CppAD::vector<size_t>, CppAD::vector<double> >(
            adPack[D].pattern_hess);
    adPack[D].match_hess = Rcpp::as<std::vector<size_t>>(matchHere);
    adPack[D].p_hess = Rcpp::as<std::vector<size_t>>(pHere);

  }

}

inline Rcpp::List sparsity(
	std::vector<GroupPack>& adPack,
	const CppAD::vector<double>& parameters,
	const bool verbose=false
	) {

	const size_t Nfun = adPack.size();
	const size_t Nparams = parameters.size();

	if(verbose) {
		Rcpp::Rcout << "sparsity Nparams " << Nparams << " Nfun " << Nfun << "\n";
	}

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

    	adPack[D].fun.for_jac_sparsity(
    		pattern_in,
    		transpose,
    		dependency,
    		internal_bool,
    		pattern_out_grad);


    	adPack[D].fun.for_hes_sparsity(
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

#endif