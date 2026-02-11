#include <Rcpp.h>

#include "adlaplace/runtime/interfaces.hpp"
#include "adlaplace/api/callables.hpp"

//' @export
// [[Rcpp::export]]
Rcpp::List getAdFun_r(
  Rcpp::List data,
  Rcpp::List config,
  const bool inner = false) {
  (void)inner;
  return getAdFun_h(data, config);
}

//' @export
// [[Rcpp::export]]
double jointLogDensNoAdfun(
  Rcpp::NumericVector parameters,
  Rcpp::List data,
  Rcpp::List config) {
  Rcpp::List adFun = getAdFun_h(data, config);
  SEXP out = adlaplace_joint_log_dens_callable()(parameters, adFun, R_NilValue);
  return Rcpp::as<double>(out);
}

//' @export
// [[Rcpp::export]]
double jointLogDens(
  Rcpp::NumericVector parameters,
  SEXP adFun,
  Rcpp::List config) {
  (void)config;
  SEXP out = adlaplace_joint_log_dens_callable()(parameters, adFun, R_NilValue);
  return Rcpp::as<double>(out);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
  Rcpp::NumericVector parameters,
  SEXP adFun,
  Rcpp::List config) {
  (void)config;
  SEXP out = adlaplace_grad_callable()(parameters, adFun, Rcpp::wrap(false), R_NilValue);
  return Rcpp::as<Rcpp::NumericVector>(out);
}

//' @export
// [[Rcpp::export]]
Rcpp::S4 hessian(
  Rcpp::NumericVector parameters,
  SEXP adFun,
  Rcpp::List config) {
  (void)config;
  SEXP out = adlaplace_hess_callable()(parameters, adFun, Rcpp::wrap(false), R_NilValue);
  return Rcpp::as<Rcpp::S4>(out);
}

//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
  Rcpp::NumericVector parameters,
  Rcpp::List data,
  Rcpp::List control,
  Rcpp::List config,
  SEXP adFun = R_NilValue) {
  SEXP adFun_local = adFun;
  if (adFun_local == R_NilValue) {
    adFun_local = getAdFun_h(data, config);
  }
  Rcpp::NumericVector gamma = config["gamma"];
  Rcpp::Environment ns = Rcpp::Environment::namespace_env("adlaplace");
  Rcpp::Function inner_opt_fn = ns["inner_opt"];
  return Rcpp::as<Rcpp::List>(inner_opt_fn(parameters, gamma, config, control, adFun_local));
}

//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt_adpack(
  Rcpp::NumericVector parameters,
  SEXP adFun,
  const Rcpp::List control,
  const Rcpp::List config,
  SEXP adFunFull = R_NilValue) {
  (void)adFunFull;
  Rcpp::NumericVector gamma = config["gamma"];
  Rcpp::Environment ns = Rcpp::Environment::namespace_env("adlaplace");
  Rcpp::Function inner_opt_fn = ns["inner_opt"];
  return Rcpp::as<Rcpp::List>(inner_opt_fn(parameters, gamma, config, control, adFun));
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
  const Rcpp::NumericVector parameters,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  const Rcpp::List config,
  SEXP adFun = R_NilValue) {
  (void)config;
  if (adFun == R_NilValue) {
    Rcpp::stop("traceHinvT requires non-NULL adFun");
  }
  SEXP out = adlaplace_trace_hinv_t_callable()(
    parameters,
    adFun,
    LinvPt,
    LinvPtColumns,
    R_NilValue
  );
  return Rcpp::as<Rcpp::NumericVector>(out);
}

//' @export
// [[Rcpp::export]]
Rcpp::List sparsity(
   Rcpp::List data,
   Rcpp::List config) {
  Rcpp::List adFun = getAdFun_h(data, config);
  return Rcpp::as<Rcpp::List>(adFun["sparsity"]);
}
