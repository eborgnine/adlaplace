#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "adlaplace/runtime/interfaces.hpp"

//' Build raw AD pack external pointer (adlaplaceExample)
//'
//' @param data Model data list.
//' @param config Model configuration list.
//' @return External pointer of class \code{adlaplace_handle_ptr}.
//' @keywords internal
// [[Rcpp::export]]
SEXP get_ad_fun_raw_example(
  Rcpp::List data,
  Rcpp::List config)
{
  return get_ad_fun_raw_h(data, config);
}
