#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "adlaplace/runtime/interfaces.hpp"

//' Backend AD builder for adlaplaceExample
//'
//' Builds model-specific AD tapes/metadata for this backend package.
//' Evaluation should be done through \pkg{adlaplace} functions
//' (\code{jointLogDens()}, \code{grad()}, \code{hess()}, \code{traceHinvT()}).
//'
//' @param data Model data list.
//' @param config Model configuration list.
//'
//' @return List with backend handle/sparsity metadata used by \pkg{adlaplace}.
//'
//' @name adlaplace_backend_cpp

//' @rdname adlaplace_backend_cpp
//' @export
// [[Rcpp::export]]
Rcpp::List getAdFun_r(
  Rcpp::List data,
  Rcpp::List config)
{
  return getAdFun_h(data, config);
}
