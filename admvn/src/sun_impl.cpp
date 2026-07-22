#include <Rcpp.h>
#include <cppad/cppad.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sun_holder.hpp"
#include "sun_tape.hpp"

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

Rcpp::List sun_result_to_list(const admvn::SunResult& res, int deriv) {
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

admvn::SunTapeBundle make_sun_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads) {
  return admvn::create_sun_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads);
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
  int n_threads = 1) {

  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);

  if (x_rows.empty()) {
    admvn::SunResult out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSunNPar, 0.0);
    out.hessian.assign(admvn::kSunNPar, std::vector<double>(admvn::kSunNPar, 0.0));
    return sun_result_to_list(out, deriv);
  }

  admvn::SunTapeBundle bundle = make_sun_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads);
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
  int n_threads = 1) {

  auto* holder = new admvn::SunTapeHolder(
    make_sun_bundle(
      as_matrix_rows(x),
      as_vector(par_seed),
      n_points,
      n_shifts,
      seed,
      n_threads),
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
  return omp_get_max_threads();
#else
  return 1;
#endif
}
