
#ifndef TRUSTOPTIMCONTROL_HPP
#define TRUSTOPTIMCONTROL_HPP

#include<Rcpp.h>

// Small helpers with defaults
inline double get_double_ctrl(const Rcpp::List& ctl,
  const char* key,
  double def) {
  return ctl.containsElementNamed(key) ? Rcpp::as<double>(ctl[key]) : def;
}

inline int get_int_ctrl(const Rcpp::List& ctl,
  const char* key,
  int def) {
  return ctl.containsElementNamed(key) ? Rcpp::as<int>(ctl[key]) : def;
}



// Bundle of trust-region control parameters
struct TrustControl {
  // Scalars
  double rad;
  double min_rad;
  double tol;
  double prec;

  int    report_freq;
  int    report_level;
  int    header_freq;
  int    report_precision;
  int    maxit;

  double contract_factor;
  double expand_factor;
  double contract_threshold;
  double expand_threshold_rad;
  double expand_threshold_ap;
  double function_scale_factor;

  int    precond_refresh_freq;
  int    precond_ID;
//  int    quasi_newton_method;  // 1 = SR1, 2 = BFGS
  int    trust_iter;

  // Construct from an Rcpp::List with defaults
  explicit TrustControl(const Rcpp::List& control)
  : rad(                  get_double_ctrl(control, "step.size",              1.0))
  , min_rad(              get_double_ctrl(control, "min.step.size",          1e-8))
  , tol(                  get_double_ctrl(control, "cg.tol",                 1e-4))
  , prec(                 get_double_ctrl(control, "grad.tol",               1e-6))
  , report_freq(          get_int_ctrl   (control, "report.freq",            0))
  , report_level(         get_int_ctrl   (control, "report.level",           0))
  , header_freq(          get_int_ctrl   (control, "header.freq",            10))
  , report_precision(     get_int_ctrl   (control, "report.precision",       6))
  , maxit(                get_int_ctrl   (control, "maxit",                  100))
  , contract_factor(      get_double_ctrl(control, "contract.factor",        0.5))
  , expand_factor(        get_double_ctrl(control, "expand.factor",          2.0))
  , contract_threshold(   get_double_ctrl(control, "contract.threshold",     0.25))
  , expand_threshold_rad( get_double_ctrl(control, "expand.threshold.rad",   0.8))
  , expand_threshold_ap(  get_double_ctrl(control, "expand.threshold.ap",    0.75))
  , function_scale_factor(get_double_ctrl(control, "function.scale.factor",  1.0))
  , precond_refresh_freq( get_int_ctrl   (control, "precond.refresh",        5))
  , precond_ID(           get_int_ctrl   (control, "precond.ID",             0))
//    , quasi_newton_method(  get_int_ctrl   (control, "quasi.newton.method",    1))  // 1 = SR1
  , trust_iter(           get_int_ctrl   (control, "trust.iter",             50))
  {}
};

#endif
