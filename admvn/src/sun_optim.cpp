#include "ompad.hpp"
#include "sun_holder.hpp"
#include "sun_tape.hpp"

#include <Rcpp.h>
#include <RcppEigen.h>

#include <CG-base.h>
#include <CG-quasi.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

inline double get_double_ctrl(const Rcpp::List& ctl, const char* key, double def) {
  if (ctl.containsElementNamed(key) && !Rf_isNull(ctl[key])) {
    return Rcpp::as<double>(ctl[key]);
  }
  return def;
}

inline int get_int_ctrl(const Rcpp::List& ctl, const char* key, int def) {
  if (ctl.containsElementNamed(key) && !Rf_isNull(ctl[key])) {
    return Rcpp::as<int>(ctl[key]);
  }
  return def;
}

struct TrustControl {
  double rad;
  double min_rad;
  double tol;
  double prec;
  int report_freq;
  int report_level;
  int header_freq;
  int report_precision;
  int maxit;
  double contract_factor;
  double expand_factor;
  double contract_threshold;
  double expand_threshold_rad;
  double expand_threshold_ap;
  double function_scale_factor;
  int precond_refresh_freq;
  int precond_ID;
  int quasi_newton_method;
  int trust_iter;

  explicit TrustControl(const Rcpp::List& control)
    : rad(get_double_ctrl(control, "step.size", 1.0)),
      min_rad(get_double_ctrl(control, "min.step.size", 1e-8)),
      tol(get_double_ctrl(control, "cg.tol", 1e-4)),
      prec(get_double_ctrl(control, "grad.tol", 1e-6)),
      report_freq(get_int_ctrl(control, "report.freq", 0)),
      report_level(get_int_ctrl(control, "report.level", 0)),
      header_freq(get_int_ctrl(control, "header.freq", 10)),
      report_precision(get_int_ctrl(control, "report.precision", 6)),
      maxit(get_int_ctrl(control, "maxit", 100)),
      contract_factor(get_double_ctrl(control, "contract.factor", 0.5)),
      expand_factor(get_double_ctrl(control, "expand.factor", 2.0)),
      contract_threshold(get_double_ctrl(control, "contract.threshold", 0.25)),
      expand_threshold_rad(get_double_ctrl(control, "expand.threshold.rad", 0.8)),
      expand_threshold_ap(get_double_ctrl(control, "expand.threshold.ap", 0.75)),
      function_scale_factor(get_double_ctrl(control, "function.scale.factor", 1.0)),
      precond_refresh_freq(get_int_ctrl(control, "precond.refresh", 5)),
      precond_ID(get_int_ctrl(control, "precond.ID", 0)),
      quasi_newton_method(get_int_ctrl(control, "quasi.newton.method", 1)),
      trust_iter(get_int_ctrl(control, "trust.iter", 50)) {}
};

// Indices of nu1, nu2, nu3 in the 21-parameter SUN vector.
constexpr int kNuBegin = 3;
constexpr int kNuEnd = 6;

std::vector<double> unconstrained_to_natural(const Eigen::VectorXd& u) {
  std::vector<double> par(static_cast<std::size_t>(u.size()));
  for (Eigen::Index i = 0; i < u.size(); ++i) {
    par[static_cast<std::size_t>(i)] = u[i];
  }
  for (int i = kNuBegin; i < kNuEnd; ++i) {
    par[static_cast<std::size_t>(i)] = std::exp(u[i]);
  }
  return par;
}

Eigen::VectorXd natural_to_unconstrained(const std::vector<double>& par) {
  Eigen::VectorXd u(static_cast<Eigen::Index>(par.size()));
  for (std::size_t i = 0; i < par.size(); ++i) {
    u[static_cast<Eigen::Index>(i)] = par[i];
  }
  for (int i = kNuBegin; i < kNuEnd; ++i) {
    const double nu = par[static_cast<std::size_t>(i)];
    if (!(nu > 0.0)) {
      Rcpp::stop("nu parameters must be positive for log-reparameterization");
    }
    u[i] = std::log(nu);
  }
  return u;
}

// Minimize -loglik(par(u)) where u has log(nu).
struct SunMleFunc {
  admvn::SunTapeBundle& bundle;
  int n_threads;

  int get_nvars() const {
    return static_cast<int>(admvn::kSunNPar);
  }

  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX>& x, double& f) {
    Eigen::VectorXd g_dummy(x.size());
    get_fdf(x, f, g_dummy);
  }

  template <class DerivedX, class DerivedG>
  void get_df(
    const Eigen::MatrixBase<DerivedX>& x,
    Eigen::MatrixBase<DerivedG>& g) {
    double f_dummy = 0.0;
    get_fdf(x, f_dummy, g);
  }

  template <class DerivedX, class DerivedG>
  void get_fdf(
    const Eigen::MatrixBase<DerivedX>& x,
    double& f,
    Eigen::MatrixBase<DerivedG>& g) {

    const std::vector<double> par = unconstrained_to_natural(x.derived());
    admvn::SunResult res = admvn::eval_sun_bundle(
      bundle, par, /*log_scale=*/true, /*deriv=*/1, n_threads,
      /*manage_parallel_scope=*/false);

    // Minimize negative log-likelihood.
    f = -res.value;
    auto& gout = g.derived();
    gout.resize(static_cast<Eigen::Index>(admvn::kSunNPar));
    for (std::size_t i = 0; i < admvn::kSunNPar; ++i) {
      gout[static_cast<Eigen::Index>(i)] = -res.gradient[i];
    }
    // Chain rule: d/d log(nu) = nu * d/d nu
    for (int i = kNuBegin; i < kNuEnd; ++i) {
      gout[i] *= par[static_cast<std::size_t>(i)];
    }
  }
};

}  // namespace

// [[Rcpp::export]]
Rcpp::List sun_mle_cpp(
  SEXP ptr,
  Rcpp::NumericVector start,
  Rcpp::List control,
  int n_threads = 0) {

  admvn::SunTapeHolder* holder = admvn::sun_holder_from_sexp(ptr);
  if (static_cast<std::size_t>(start.size()) != admvn::kSunNPar) {
    Rcpp::stop("start must have length 21");
  }

  std::vector<double> start_nat = Rcpp::as<std::vector<double>>(start);
  Eigen::VectorXd start_u = natural_to_unconstrained(start_nat);

  int threads = n_threads;
  if (threads <= 0) {
    threads = holder->bundle.n_threads;
  }
  if (threads < 1) {
    threads = 1;
  }

  TrustControl ctl(control);
  SunMleFunc fun{holder->bundle, threads};

  using Tvec = Eigen::VectorXd;
  using THess = Eigen::MatrixXd;
  using TPreLLt = Eigen::LLT<Eigen::MatrixXd>;

  CppadParallelScope parallel_scope(static_cast<std::size_t>(threads));

  Trust_CG_Optimizer<Tvec, SunMleFunc, THess, TPreLLt> opt(
    fun,
    start_u,
    ctl.rad,
    ctl.min_rad,
    ctl.tol,
    ctl.prec,
    ctl.report_freq,
    ctl.report_level,
    ctl.header_freq,
    ctl.report_precision,
    ctl.maxit,
    ctl.contract_factor,
    ctl.expand_factor,
    ctl.contract_threshold,
    ctl.expand_threshold_rad,
    ctl.expand_threshold_ap,
    ctl.function_scale_factor,
    ctl.precond_refresh_freq,
    ctl.precond_ID,
    ctl.quasi_newton_method,
    ctl.trust_iter);

  opt.run();

  Tvec sol(admvn::kSunNPar);
  Tvec grad_u(admvn::kSunNPar);
  double fval = NA_REAL;
  double radius = NA_REAL;
  int iterations = NA_INTEGER;
  MB_Status status = opt.get_current_state(sol, fval, grad_u, iterations, radius);

  std::vector<double> par_hat = unconstrained_to_natural(sol);
  // Report positive log-likelihood and natural-scale gradient.
  admvn::SunResult at_hat = admvn::eval_sun_bundle(
    holder->bundle, par_hat, true, 1, threads, /*manage_parallel_scope=*/false);

  Rcpp::NumericVector par_out(par_hat.begin(), par_hat.end());
  if (!Rf_isNull(start.names())) {
    par_out.attr("names") = start.names();
  }

  return Rcpp::List::create(
    Rcpp::Named("par") = par_out,
    Rcpp::Named("value") = at_hat.value,
    Rcpp::Named("gradient") = Rcpp::NumericVector(
      at_hat.gradient.begin(), at_hat.gradient.end()),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("radius") = radius,
    Rcpp::Named("status") = static_cast<int>(status),
    Rcpp::Named("f_min") = fval,
    Rcpp::Named("n_threads") = threads);
}
