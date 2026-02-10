#ifndef ADLAPLACE_CREATORS_CALLABLES_HPP
#define ADLAPLACE_CREATORS_CALLABLES_HPP

#include <Rcpp.h>
#include <R_ext/Rdynload.h>

using adlaplace_gh_callable_t = SEXP (*)(SEXP, SEXP, SEXP, SEXP);
using adlaplace_joint_callable_t = SEXP (*)(SEXP, SEXP, SEXP);
using adlaplace_trace_hinv_t_callable_t = SEXP (*)(SEXP, SEXP, SEXP, SEXP, SEXP);

inline DL_FUNC lookup_adlaplace_callable(const char* symbol) {
  // Ensure adlaplace has registered its C-callables before lookup.
  Rcpp::Environment adlaplace_ns = Rcpp::Environment::namespace_env("adlaplace");
  if (adlaplace_ns.exists("register_callables")) {
    Rcpp::Function reg = adlaplace_ns["register_callables"];
    reg();
  }
  DL_FUNC fn = R_GetCCallable("adlaplace", symbol);
  if (fn == nullptr) {
    Rcpp::stop("Failed to load C-callable '%s' from package 'adlaplace'", symbol);
  }
  return fn;
}

inline adlaplace_gh_callable_t adlaplace_grad_callable() {
  static adlaplace_gh_callable_t fn = nullptr;
  if (fn == nullptr) {
    fn = reinterpret_cast<adlaplace_gh_callable_t>(
      lookup_adlaplace_callable("adlaplace_grad_c")
    );
  }
  return fn;
}

inline adlaplace_gh_callable_t adlaplace_hess_callable() {
  static adlaplace_gh_callable_t fn = nullptr;
  if (fn == nullptr) {
    fn = reinterpret_cast<adlaplace_gh_callable_t>(
      lookup_adlaplace_callable("adlaplace_hess_c")
    );
  }
  return fn;
}

inline adlaplace_joint_callable_t adlaplace_joint_log_dens_callable() {
  static adlaplace_joint_callable_t fn = nullptr;
  if (fn == nullptr) {
    fn = reinterpret_cast<adlaplace_joint_callable_t>(
      lookup_adlaplace_callable("adlaplace_jointLogDens_c")
    );
  }
  return fn;
}

inline adlaplace_trace_hinv_t_callable_t adlaplace_trace_hinv_t_callable() {
  static adlaplace_trace_hinv_t_callable_t fn = nullptr;
  if (fn == nullptr) {
    fn = reinterpret_cast<adlaplace_trace_hinv_t_callable_t>(
      lookup_adlaplace_callable("adlaplace_traceHinvT_c")
    );
  }
  return fn;
}

#endif
