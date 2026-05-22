#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "adlaplace/api/register.hpp"

//' Build raw AD pack external pointer (hpolcc)
//'
//' @param data Model data list.
//' @param config Model configuration list.
//' @return External pointer of class \code{adlaplace_handle_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_hpolcc(
  Rcpp::List data,
  Rcpp::List config) {

  return get_ad_fun_raw_h(data, config);
}
