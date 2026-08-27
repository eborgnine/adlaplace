#include <Rcpp.h>
#include <cppad/cppad.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sun_holder.hpp"
#include "sun_tape.hpp"
#include "sun22_holder.hpp"
#include "sun22_tape.hpp"
#include "sun32_holder.hpp"
#include "sun32_tape.hpp"
#include "sun42_holder.hpp"
#include "sun42_tape.hpp"
#include "sun43_holder.hpp"
#include "sun43_tape.hpp"
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

Rcpp::List sun22_result_to_list(const admvn::Sun22Result& res, int deriv) {
  return sun_like_result_to_list(res, deriv);
}

Rcpp::List sun32_result_to_list(const admvn::Sun32Result& res, int deriv) {
  return sun_like_result_to_list(res, deriv);
}

Rcpp::List sun42_result_to_list(const admvn::Sun42Result& res, int deriv) {
  return sun_like_result_to_list(res, deriv);
}

Rcpp::List sun43_result_to_list(const admvn::Sun43Result& res, int deriv) {
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
  const std::vector<double>& weights = {},
  admvn::SunParMap par_map = admvn::SunParMap::kBlockChol) {
  return admvn::create_sun_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights,
    par_map);
}

admvn::Sun22TapeBundle make_sun22_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {},
  admvn::Sun22ParMap par_map = admvn::Sun22ParMap::kBlockChol) {
  return admvn::create_sun22_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights,
    par_map);
}

admvn::Sun32TapeBundle make_sun32_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {},
  admvn::Sun32ParMap par_map = admvn::Sun32ParMap::kHyperspherical) {
  return admvn::create_sun32_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights,
    par_map);
}

admvn::Sun42TapeBundle make_sun42_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {},
  admvn::Sun42ParMap par_map = admvn::Sun42ParMap::kHyperspherical) {
  return admvn::create_sun42_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights,
    par_map);
}

admvn::Sun43TapeBundle make_sun43_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {},
  admvn::Sun43ParMap par_map = admvn::Sun43ParMap::kHyperspherical) {
  return admvn::create_sun43_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights,
    par_map);
}

admvn::Sun44TapeBundle make_sun44_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  int n_points,
  int n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights = {},
  admvn::Sun44ParMap par_map = admvn::Sun44ParMap::kBlockChol) {
  return admvn::create_sun44_bundle(
    x_rows,
    par_seed,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed,
    n_threads,
    weights,
    par_map);
}

void require_sun22_dims(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& par) {
  if (x.ncol() != static_cast<int>(admvn::kSun22D)) {
    Rcpp::stop("x must have ncol == 2 for SUN(2,2)");
  }
  if (par.size() != static_cast<int>(admvn::kSun22NPar)) {
    Rcpp::stop("par must have length 10 for SUN(2,2)");
  }
}

void require_sun32_dims(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& par) {
  if (x.ncol() != static_cast<int>(admvn::kSun32D)) {
    Rcpp::stop("x must have ncol == 3 for SUN(3,2)");
  }
  if (par.size() != static_cast<int>(admvn::kSun32NPar)) {
    Rcpp::stop("par must have length 16 for SUN(3,2)");
  }
}

void require_sun42_dims(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& par) {
  if (x.ncol() != static_cast<int>(admvn::kSun42D)) {
    Rcpp::stop("x must have ncol == 4 for SUN(4,2)");
  }
  if (par.size() != static_cast<int>(admvn::kSun42NPar)) {
    Rcpp::stop("par must have length 23 for SUN(4,2)");
  }
}

void require_sun43_dims(const Rcpp::NumericMatrix& x, const Rcpp::NumericVector& par) {
  if (x.ncol() != static_cast<int>(admvn::kSun43D)) {
    Rcpp::stop("x must have ncol == 4 for SUN(4,3)");
  }
  if (par.size() != static_cast<int>(admvn::kSun43NPar)) {
    Rcpp::stop("par must have length 29 for SUN(4,3)");
  }
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
Rcpp::List dsun_hs_cpp(
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
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v,
    admvn::SunParMap::kHyperspherical);
  admvn::SunResult res = admvn::eval_sun_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun_hs_fun_create_cpp(
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
      w_v,
      admvn::SunParMap::kHyperspherical),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun_tape_ptr"));
  return out;
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
Rcpp::List dsun22_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun22_dims(x, par);

  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }

  if (x_rows.empty()) {
    admvn::Sun22Result out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSun22NPar, 0.0);
    out.hessian.assign(admvn::kSun22NPar, std::vector<double>(admvn::kSun22NPar, 0.0));
    return sun22_result_to_list(out, deriv);
  }

  admvn::Sun22TapeBundle bundle = make_sun22_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v);
  admvn::Sun22Result res = admvn::eval_sun22_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun22_result_to_list(res, deriv);
}

// [[Rcpp::export]]
Rcpp::List dsun22_hs_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun22_dims(x, par);
  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  if (x_rows.empty()) {
    admvn::Sun22Result out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSun22NPar, 0.0);
    out.hessian.assign(
      admvn::kSun22NPar, std::vector<double>(admvn::kSun22NPar, 0.0));
    return sun22_result_to_list(out, deriv);
  }
  admvn::Sun22TapeBundle bundle = make_sun22_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v,
    admvn::Sun22ParMap::kHyperspherical);
  admvn::Sun22Result res = admvn::eval_sun22_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun22_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun22_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun22_dims(x, par_seed);

  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }

  auto* holder = new admvn::Sun22TapeHolder(
    make_sun22_bundle(
      as_matrix_rows(x),
      as_vector(par_seed),
      n_points,
      n_shifts,
      seed,
      n_threads,
      w_v),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun22_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun22_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
SEXP dsun22_hs_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun22_dims(x, par_seed);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  auto* holder = new admvn::Sun22TapeHolder(
    make_sun22_bundle(
      as_matrix_rows(x), as_vector(par_seed), n_points, n_shifts, seed,
      n_threads, w_v, admvn::Sun22ParMap::kHyperspherical),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun22_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun22_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List dsun22_fun_eval_cpp(
  SEXP ptr,
  Rcpp::Nullable<Rcpp::NumericVector> par = R_NilValue,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0) {

  auto* holder = admvn::sun22_holder_from_sexp(ptr);

  std::vector<double> par_v;
  if (par.isNotNull()) {
    par_v = as_vector(par.get());
  } else {
    par_v = holder->seed_par;
  }
  if (par_v.size() != admvn::kSun22NPar) {
    Rcpp::stop("par must have length 10 for SUN(2,2)");
  }

  admvn::Sun22Result res = admvn::eval_sun22_bundle(
    holder->bundle, par_v, log_scale, deriv, n_threads);
  return sun22_result_to_list(res, deriv);
}

// [[Rcpp::export]]
Rcpp::List dsun32_hs_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun32_dims(x, par);
  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  if (x_rows.empty()) {
    admvn::Sun32Result out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSun32NPar, 0.0);
    out.hessian.assign(
      admvn::kSun32NPar, std::vector<double>(admvn::kSun32NPar, 0.0));
    return sun32_result_to_list(out, deriv);
  }
  admvn::Sun32TapeBundle bundle = make_sun32_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v,
    admvn::Sun32ParMap::kHyperspherical);
  admvn::Sun32Result res = admvn::eval_sun32_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun32_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun32_hs_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun32_dims(x, par_seed);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  auto* holder = new admvn::Sun32TapeHolder(
    make_sun32_bundle(
      as_matrix_rows(x), as_vector(par_seed), n_points, n_shifts, seed,
      n_threads, w_v, admvn::Sun32ParMap::kHyperspherical),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun32_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun32_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List dsun32_fun_eval_cpp(
  SEXP ptr,
  Rcpp::Nullable<Rcpp::NumericVector> par = R_NilValue,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0) {

  auto* holder = admvn::sun32_holder_from_sexp(ptr);

  std::vector<double> par_v;
  if (par.isNotNull()) {
    par_v = as_vector(par.get());
  } else {
    par_v = holder->seed_par;
  }
  if (par_v.size() != admvn::kSun32NPar) {
    Rcpp::stop("par must have length 16 for SUN(3,2)");
  }

  admvn::Sun32Result res = admvn::eval_sun32_bundle(
    holder->bundle, par_v, log_scale, deriv, n_threads);
  return sun32_result_to_list(res, deriv);
}

// [[Rcpp::export]]
Rcpp::List dsun42_hs_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun42_dims(x, par);
  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  if (x_rows.empty()) {
    admvn::Sun42Result out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSun42NPar, 0.0);
    out.hessian.assign(
      admvn::kSun42NPar, std::vector<double>(admvn::kSun42NPar, 0.0));
    return sun42_result_to_list(out, deriv);
  }
  admvn::Sun42TapeBundle bundle = make_sun42_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v,
    admvn::Sun42ParMap::kHyperspherical);
  admvn::Sun42Result res = admvn::eval_sun42_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun42_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun42_hs_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun42_dims(x, par_seed);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  auto* holder = new admvn::Sun42TapeHolder(
    make_sun42_bundle(
      as_matrix_rows(x), as_vector(par_seed), n_points, n_shifts, seed,
      n_threads, w_v, admvn::Sun42ParMap::kHyperspherical),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun42_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun42_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List dsun42_fun_eval_cpp(
  SEXP ptr,
  Rcpp::Nullable<Rcpp::NumericVector> par = R_NilValue,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0) {

  auto* holder = admvn::sun42_holder_from_sexp(ptr);

  std::vector<double> par_v;
  if (par.isNotNull()) {
    par_v = as_vector(par.get());
  } else {
    par_v = holder->seed_par;
  }
  if (par_v.size() != admvn::kSun42NPar) {
    Rcpp::stop("par must have length 23 for SUN(4,2)");
  }

  admvn::Sun42Result res = admvn::eval_sun42_bundle(
    holder->bundle, par_v, log_scale, deriv, n_threads);
  return sun42_result_to_list(res, deriv);
}

// [[Rcpp::export]]
Rcpp::List dsun43_hs_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par,
  bool log_scale = true,
  int deriv = 0,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun43_dims(x, par);
  const auto x_rows = as_matrix_rows(x);
  const std::vector<double> par_v = as_vector(par);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  if (x_rows.empty()) {
    admvn::Sun43Result out;
    out.value = log_scale ? 0.0 : 1.0;
    out.gradient.assign(admvn::kSun43NPar, 0.0);
    out.hessian.assign(
      admvn::kSun43NPar, std::vector<double>(admvn::kSun43NPar, 0.0));
    return sun43_result_to_list(out, deriv);
  }
  admvn::Sun43TapeBundle bundle = make_sun43_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v,
    admvn::Sun43ParMap::kHyperspherical);
  admvn::Sun43Result res = admvn::eval_sun43_bundle(
    bundle, par_v, log_scale, deriv, n_threads);
  return sun43_result_to_list(res, deriv);
}

// [[Rcpp::export]]
SEXP dsun43_hs_fun_create_cpp(
  Rcpp::NumericMatrix x,
  Rcpp::NumericVector par_seed,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1,
  int n_threads = 1,
  Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {

  require_sun43_dims(x, par_seed);
  std::vector<double> w_v;
  if (weights.isNotNull()) {
    w_v = as_vector(weights.get());
  }
  auto* holder = new admvn::Sun43TapeHolder(
    make_sun43_bundle(
      as_matrix_rows(x), as_vector(par_seed), n_points, n_shifts, seed,
      n_threads, w_v, admvn::Sun43ParMap::kHyperspherical),
    as_vector(par_seed));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, admvn::sun43_holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_sun43_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List dsun43_fun_eval_cpp(
  SEXP ptr,
  Rcpp::Nullable<Rcpp::NumericVector> par = R_NilValue,
  bool log_scale = true,
  int deriv = 0,
  int n_threads = 0) {

  auto* holder = admvn::sun43_holder_from_sexp(ptr);

  std::vector<double> par_v;
  if (par.isNotNull()) {
    par_v = as_vector(par.get());
  } else {
    par_v = holder->seed_par;
  }
  if (par_v.size() != admvn::kSun43NPar) {
    Rcpp::stop("par must have length 29 for SUN(4,3)");
  }

  admvn::Sun43Result res = admvn::eval_sun43_bundle(
    holder->bundle, par_v, log_scale, deriv, n_threads);
  return sun43_result_to_list(res, deriv);
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
Rcpp::List dsun44_hs_cpp(
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
    out.hessian.assign(
      admvn::kSun44NPar, std::vector<double>(admvn::kSun44NPar, 0.0));
    return sun44_result_to_list(out, deriv);
  }
  admvn::Sun44TapeBundle bundle = make_sun44_bundle(
    x_rows, par_v, n_points, n_shifts, seed, n_threads, w_v,
    admvn::Sun44ParMap::kHyperspherical);
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
SEXP dsun44_hs_fun_create_cpp(
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
      as_matrix_rows(x), as_vector(par_seed), n_points, n_shifts, seed,
      n_threads, w_v, admvn::Sun44ParMap::kHyperspherical),
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
