#include <Rcpp.h>
#include <cppad/cppad.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sun_holder.hpp"
#include "sun_tape.hpp"
#include "sun44_holder.hpp"
#include "sun44_tape.hpp"

namespace {

std::vector<double> as_vector(const Rcpp::NumericVector& x) {
  return Rcpp::as<std::vector<double>>(x);
}

std::vector<std::vector<double>> as_matrix_rows(const Rcpp::NumericMatrix& m) {
  const int nr = m.nrow();
  const int nc = m.ncol();
  std::vector<std::vector<double>> out(nr, std::vector<double>(nc));
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      out[i][j] = m(i, j);
    }
  }
  return out;
}

template <typename Result>
Rcpp::List sun_like_result_to_list(const Result& res, int deriv) {
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("value") = res.value,
    Rcpp::Named("error") = res.error
  );
  if (deriv >= 1) {
    out["gradient"] = Rcpp::NumericVector(res.gradient.begin(), res.gradient.end());
  }
  if (deriv >= 2) {
    const int p = static_cast<int>(res.hessian.size());
    Rcpp::NumericMatrix H(p, p);
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j < p; ++j) {
        H(i, j) = res.hessian[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
      }
    }
    out["hessian"] = H;
  }
  return out;
}

Rcpp::List sun_result_to_list(const admvn::SunResult& res, int deriv) {
  return sun_like_result_to_list(res, deriv);
}

Rcpp::List sun44_result_to_list(const admvn::Sun44Result& res, int deriv) {
  return sun_like_result_to_list(res, deriv);
}

admvn::SunTapeBundle make_sun_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {}) {
  return admvn::create_sun_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights);
}

admvn::Sun44TapeBundle make_sun44_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {}) {
  return admvn::create_sun44_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights);
}

void require_sun44_dims(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& par) {
  if (x.ncol() != static_cast<int>(admvn::kSun44D)) {
    Rcpp::stop("x must have ncol == 4 for SUN(4,4)");
  }
  if (par.size() != static_cast<int>(admvn::kSun44NPar)) {
    Rcpp::stop("par must have length 36 for SUN(4,4)");
  }
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List dsun_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }

  if (x_rows.empty()) {
    admvn::SunResult out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSunNPar, 0.0);
    out.hessian.assign(admvn::kSunNPar, std::vector<double>(admvn::kSunNPar, 0.0));
    return sun_result_to_list(out, deriv);
  }

  admvn::SunTapeBundle bundle = make_sun_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v);
  admvn::SunResult res = admvn::eval_sun_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }

  auto* holder = new admvn::SunTapeHolder(
    make_sun_bundle(
      as_matrix_rows(x),
      as_vector(par_seed),
      n_points,
      n_shifts,
      seed,
      n_threads,
      w_v),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List dsun_fun_eval_cpp(
  SEXP ptr,
  Rcpp::Nullable<Rcpp::NumericVector> par = R_NilValue,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0) {

  auto* holder = admvn::sun_holder_from_sexp(ptr);

  std::vector<double> par_v;
  if (par.isNotNull()) {
    par_v = as_vector(par.get());
  } else {
    par_v = holder->seed_par;
  }

  admvn::SunResult res = admvn::eval_sun_bundle(
    holder->bundle, par_v, log_scale, deriv, n_threads);
  return sun_result_to_list(res, deriv);
}

// [[Rcpp::export]]
int dsun_n_threads_default_cpp() {
#ifdef _OPENMP
  // Prefer hardware concurrency over omp_get_max_threads(), which tracks the
  // last omp_set_num_threads() (often 1 after CppAD parallel teardown).
  const int n = omp_get_num_procs();
  return n > 0 ? n : 1;
#else
  return 1;
#endif
}

// [[Rcpp::export]]
Rcpp::List dsun44_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun44_dims(x, par);

  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }

  if (x_rows.empty()) {
    admvn::Sun44Result out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSun44NPar, 0.0);
    out.hessian.assign(admvn::kSun44NPar, std::vector<double>(admvn::kSun44NPar, 0.0));
    return sun44_result_to_list(out, deriv);
  }

  admvn::Sun44TapeBundle bundle = make_sun44_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v);
  admvn::Sun44Result res = admvn::eval_sun44_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun44_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun44_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun44_dims(x, par_seed);

  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }

  auto* holder = new admvn::Sun44TapeHolder(
    make_sun44_bundle(
      as_matrix_rows(x),
      as_vector(par_seed),
      n_points,
      n_shifts,
      seed,
      n_threads,
      w_v),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun44_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun44_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List dsun44_fun_eval_cpp(
  SEXP ptr,
  Rcpp::Nullable<Rcpp::NumericVector> par = R_NilValue,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0) {

  auto* holder = admvn::sun44_holder_from_sexp(ptr);

  std::vector<double> par_v;
  if (par.isNotNull()) {
    par_v = as_vector(par.get());
  } else {
    par_v = holder->seed_par;
  }
  if (par_v.size() != admvn::kSun44NPar) {
    Rcpp::stop("par must have length 36 for SUN(4,4)");
  }

  admvn::Sun44Result res = admvn::eval_sun44_bundle(
    holder->bundle, par_v, log_scale, deriv, n_threads);
  return sun44_result_to_list(res, deriv);
}
