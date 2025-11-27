#ifndef ADFUN_HPP
#define ADFUN_HPP

#include "adlaplace/adlaplace.hpp"

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
  const Rcpp::List& sparsity)
{

  const size_t Nfun = adPack.size();
  const size_t Nparams = adPack[0].fun.Domain();

  if(sparsity.size() != Nfun) {
    Rcpp::warning("sparsity size smaller than adPack ", sparsity.size(), " ", Nfun);
  }
  for(size_t D=0;D<Nfun;D++) {
    const Rcpp::List sparseHere = sparsity[D];
    const Rcpp::IntegerVector sparseGrad = sparseHere["grad"];
    const Rcpp::IntegerVector rowOut = sparseHere["i"];
    const Rcpp::IntegerVector colOut = sparseHere["j"];
    const Rcpp::IntegerVector matchHere = sparseHere["match"];
    const Rcpp::IntegerVector pHere = sparseHere["p"];

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

#endif
