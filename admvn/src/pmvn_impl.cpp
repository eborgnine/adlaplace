#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include "mvn_tape.hpp"

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

class MvnTapeHolder {
public:
  admvn::MvnSetup setup;
  CppAD::ADFun<double> fun;

  explicit MvnTapeHolder(admvn::MvnSetup s)
    : setup(std::move(s)), fun(admvn::build_mvn_tape(setup)) {}
};

void holder_finalizer(SEXP ptr) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    return;
  }
  delete static_cast<MvnTapeHolder*>(R_ExternalPtrAddr(ptr));
  R_ClearExternalPtr(ptr);
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

  admvn::MvnSetup setup = admvn::prepare_mvn_setup(
    lower_v, upper_v, mean_v, sigma_v,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed);

  admvn::MvnResult res = admvn::eval_mvn_tape(setup, upper_v, mean_v, nullptr);

  Rcpp::NumericMatrix H(
    static_cast<int>(res.hessian.size()),
    res.hessian.empty() ? 0 : static_cast<int>(res.hessian[0].size()));
  for (std::size_t i = 0; i < res.hessian.size(); ++i) {
    for (std::size_t j = 0; j < res.hessian[i].size(); ++j) {
      H(static_cast<int>(i), static_cast<int>(j)) = res.hessian[i][j];
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = res.value,
    Rcpp::Named("error") = res.error,
    Rcpp::Named("gradient") = Rcpp::NumericVector(res.gradient.begin(), res.gradient.end()),
    Rcpp::Named("hessian") = H
  );
}

// [[Rcpp::export]]
SEXP pmvn_fun_create_cpp(
  Rcpp::NumericVector lower,
  Rcpp::NumericVector mean,
  Rcpp::NumericMatrix sigma,
  int n_points = 1021,
  int n_shifts = 8,
  unsigned int seed = 1) {

  const std::vector<double> lower_v = as_vector(lower);
  const std::vector<double> mean_v = as_vector(mean);
  const auto sigma_v = as_matrix(sigma);

  std::vector<double> upper_v(lower_v.size());
  for (std::size_t i = 0; i < upper_v.size(); ++i) {
    upper_v[i] = mean_v[i];
  }

  admvn::MvnSetup setup = admvn::prepare_mvn_setup(
    lower_v, upper_v, mean_v, sigma_v,
    static_cast<std::size_t>(n_points),
    static_cast<std::size_t>(n_shifts),
    seed);

  auto* holder = new MvnTapeHolder(std::move(setup));
  SEXP out = R_MakeExternalPtr(holder, R_NilValue, R_NilValue);
  R_RegisterCFinalizerEx(out, holder_finalizer, TRUE);
  Rf_setAttrib(out, R_ClassSymbol, Rf_mkString("admvn_tape_ptr"));
  return out;
}

// [[Rcpp::export]]
Rcpp::List pmvn_fun_eval_cpp(SEXP ptr, Rcpp::NumericVector upper) {
  if (R_ExternalPtrAddr(ptr) == nullptr) {
    Rcpp::stop("invalid admvn tape pointer");
  }
  auto* holder = static_cast<MvnTapeHolder*>(R_ExternalPtrAddr(ptr));
  const std::vector<double> upper_v = as_vector(upper);
  admvn::MvnResult res = admvn::eval_mvn_tape(
    holder->setup, upper_v, holder->setup.mean, &holder->fun);

  Rcpp::NumericMatrix H(
    static_cast<int>(res.hessian.size()),
    res.hessian.empty() ? 0 : static_cast<int>(res.hessian[0].size()));
  for (std::size_t i = 0; i < res.hessian.size(); ++i) {
    for (std::size_t j = 0; j < res.hessian[i].size(); ++j) {
      H(static_cast<int>(i), static_cast<int>(j)) = res.hessian[i][j];
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("value") = res.value,
    Rcpp::Named("error") = res.error,
    Rcpp::Named("gradient") = Rcpp::NumericVector(res.gradient.begin(), res.gradient.end()),
    Rcpp::Named("hessian") = H
  );
}
