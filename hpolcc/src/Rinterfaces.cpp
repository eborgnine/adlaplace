#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "adlaplace/api/register.hpp"

//' Backend AD builder for hpolcc
//'
//' Builds model-specific AD tapes/metadata for this backend package.
//' Evaluation should be done through \pkg{adlaplace} functions
//' (\code{jointLogDens()}, \code{grad()}, \code{hess()}, \code{traceHinvT()}).
//'
//' @param data Model data list.
//' @param config Model configuration list.
//'
//' @return External pointer AD handle used by \pkg{adlaplace}.
//'
//' @name adlaplace_backend_cpp

//' @rdname adlaplace_backend_cpp
//' @export
// [[Rcpp::export]]
SEXP getAdFun_r(
  Rcpp::List data,
  Rcpp::List config) {

  return get_ad_fun_raw_h(data, config);
}
