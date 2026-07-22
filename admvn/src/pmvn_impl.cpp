#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "mvn_analytic.hpp"
#include "mvn_tape.hpp"

#include <limits>

namespace {

std::vector<double> as_vector(const Rcpp::NumericVector& x) {
  return Rcpp::as<std::vector<double>>(x);
}

std::vector<std::vector<double>> as_matrix(const Rcpp::NumericMatrix& m) {
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

Rcpp::List result_to_list(const admvn::MvnResult& res, bool include_mean_grad) {
  Rcpp::NumericMatrix H(
    static_cast<int>(res.hessian.size()),
    res.hessian.empty() ? 0 : static_cast<int>(res.hessian[0].size()));
  for (std::size_t i = 0; i < res.hessian.size(); ++i) {
    for (std::size_t j = 0; j < res.hessian[i].size(); ++j) {
      H(static_cast<int>(i), static_cast<int>(j)) = res.hessian[i][j];
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("value") = res.value,
    Rcpp::Named("error") = res.error,
    Rcpp::Named("gradient") = Rcpp::NumericVector(res.gradient.begin(), res.gradient.end()),
    Rcpp::Named("hessian") = H
  );

  if (include_mean_grad) {
    out["gradient_mean"] = Rcpp::NumericVector(
      res.gradient_mean.begin(), res.gradient_mean.end());
  }
  return out;
}

class MvnTapeHolder {
public:
  admvn::MvnTape tape;
  std::vector<double> seed_mean;
  std::vector<std::vector<double>> seed_sigma;

  MvnTapeHolder(
    admvn::MvnTape t,
    std::vector<double> mean,
    std::vector<std::vector<double>> sigma)
    : tape(std::move(t)),
      seed_mean(std::move(mean)),
      seed_sigma(std::move(sigma)) {}
};

void holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<MvnTapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
}

admvn::MvnTape make_tape(
  const std::vector<double>& lower,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  int n_points,
  int n_shifts,
  unsigned int seed) {

  return admvn::create_mvn_tape(
    lower, upper, mean, sigma,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List pmvn_cpp(
  Rcpp::NumericVector upper,
  Rcpp::NumericVector lower,
  Rcpp::NumericVector mean,
  Rcpp::NumericMatrix sigma,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1) {

  const std::vector<double> upper_v = as_vector(upper);
  const std::vector<double> lower_v = as_vector(lower);
  const std::vector<double> mean_v = as_vector(mean);
  const auto sigma_v = as_matrix(sigma);

  admvn::MvnTape tape = make_tape(
    lower_v, upper_v, mean_v, sigma_v, n_points, n_shifts, seed);
  admvn::MvnResult res = admvn::eval_mvn_tape(tape, upper_v, mean_v, sigma_v, true);
  return result_to_list(res, false);
}

// [[Rcpp::export]]
Rcpp::List pack_genz_ch_cpp(
  Rcpp::NumericMatrix sigma,
  Rcpp::IntegerVector perm) {

  const auto sigma_v = as_matrix(sigma);
  std::vector<int> perm_v(perm.begin(), perm.end());
  const admvn::GenzPack genz = admvn::pack_genz_ch(sigma_v, perm_v);

  const int n = static_cast<int>(genz.scale.size());
  Rcpp::NumericMatrix ch(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ch(i, j) = (j <= i) ? genz.ch[i][j] : 0.0;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("scale") = Rcpp::NumericVector(genz.scale.begin(), genz.scale.end()),
    Rcpp::Named("ch") = ch
  );
}

// [[Rcpp::export]]
SEXP pmvn_fun_create_cpp(
  Rcpp::NumericVector lower,
  Rcpp::NumericVector mean,
  Rcpp::NumericMatrix sigma,
  Rcpp::Nullable<Rcpp::NumericVector> upper_seed = R_NilValue,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1) {

  const std::vector<double> lower_v = as_vector(lower);
  const std::vector<double> mean_v = as_vector(mean);
  const auto sigma_v = as_matrix(sigma);

  std::vector<double> upper_v(mean_v.size());
  if (upper_seed.isNotNull()) {
    upper_v = as_vector(upper_seed.get());
  } else {
    upper_v = mean_v;
  }

  auto* holder = new MvnTapeHolder(
    make_tape(lower_v, upper_v, mean_v, sigma_v, n_points, n_shifts, seed),
    mean_v,
    sigma_v);
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::IntegerVector pmvn_fun_perm_cpp(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn tape pointer");
  }
  auto* holder = static_cast<MvnTapeHolder*>(R_ExternalPtrAddr(ptr));
  return Rcpp::IntegerVector(
    holder->tape.perm.begin(), holder->tape.perm.end());
}

// [[Rcpp::export]]
Rcpp::List pmvn_fun_eval_cpp(
  SEXP ptr,
  Rcpp::NumericVector upper,
  Rcpp::Nullable<Rcpp::NumericVector> mean = R_NilValue,
  Rcpp::Nullable<Rcpp::NumericMatrix> sigma = R_NilValue,
  bool inner = true) {

  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn tape pointer");
  }
  auto* holder = static_cast<MvnTapeHolder*>(R_ExternalPtrAddr(ptr));
  const std::vector<double> upper_v = as_vector(upper);

  std::vector<double> mean_v;
  if (mean.isNotNull()) {
    mean_v = as_vector(mean.get());
  } else {
    mean_v = holder->seed_mean;
  }

  std::vector<std::vector<double>> sigma_v;
  if (sigma.isNotNull()) {
    sigma_v = as_matrix(sigma.get());
  } else {
    sigma_v = holder->seed_sigma;
  }

  admvn::MvnResult res = admvn::eval_mvn_tape(
    holder->tape, upper_v, mean_v, sigma_v, inner);
  return result_to_list(res, !inner);
}

// [[Rcpp::export]]
Rcpp::NumericVector pmvn_domain_grad_cpp(
  Rcpp::NumericVector upper,
  Rcpp::NumericVector mean,
  Rcpp::NumericVector scale,
  Rcpp::NumericMatrix ch,
  Rcpp::Nullable<Rcpp::NumericVector> lower = R_NilValue) {

  const std::size_t n = static_cast<std::size_t>(upper.size());
  admvn::GenzPack genz;
  genz.scale = as_vector(scale);
  genz.ch.assign(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      genz.ch[i][j] = ch(static_cast<int>(i), static_cast<int>(j));
    }
  }
  admvn::MvnTape tape;
  tape.n = n;
  tape.n_domain = admvn::domain_size(n);
  if (lower.isNotNull()) {
    tape.lower = as_vector(lower.get());
    if (tape.lower.size() == 1) {
      tape.lower.assign(n, tape.lower[0]);
    }
  } else {
    tape.lower.assign(n, -std::numeric_limits<double>::infinity());
  }
  tape.value_only = true;
  std::vector<double> g = admvn::eval_mvn_domain_grad_auto(
    tape, as_vector(upper), as_vector(mean), genz);
  return Rcpp::NumericVector(g.begin(), g.end());
}
