#include <Rcpp.h>
#include <Rinternals.h>

#include "adlaplace/density_data.hpp"
#include "adlaplace/ad_pack.hpp"
#include "adlaplace/ad_pack_random.hpp"
#include "adlaplace/extension.hpp"
#include "adlaplace/quadform_atomic.hpp"

#include "fem_logdet_atomic.hpp"
#include "fem_ssq_analytic.hpp"

#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

// paste0("const double TWOPI = ", format(2 * Rmpfr::Const("pi", prec =
// 120)),";")
const double PIx4 = 12.566370614359172953850573533118011539;
const double PIx8 = 25.132741228718345907701147066236023078;
// paste0("const double ONEHALFLOGTWOPI = ", format(0.5 * log(2 *
//   Rmpfr::Const("pi", prec = 120))), ";")
const double ONEHALFLOGTWOPI = 0.91893853320467274178032973640561763976;

// Declared by ADLAPLACE_DEFINE_BACKEND in backend.cpp
extern ad_shard *adlaplace_fem_make_shard(AdTape &&);

// Defined in backend.cpp, where FemSsqShard can see EvalShard. Builds the
// same ad_pack as packs_to_ad_fun but with the analytic ssq shard, which
// needs constructor data the fixed ShardFactory signature cannot carry.
extern ad_pack *fem_ssq_packs_to_ad_fun(std::vector<AdTape> &&, std::size_t,
                                        std::size_t,
                                        std::shared_ptr<const femssq::Data>);

namespace {

struct FemCholPayload {
  CscMatrix Q;   // pattern-only (x empty)
  std::vector<double> C_x;
  std::vector<double> G_x;
  std::vector<double> G2_x;
  std::vector<double> G3_x;
  std::vector<int> perm;
  CscMatrix L1;  // pattern-only from chol$L1
  int alpha = 2;
};

void fill_fem_Q_pattern(CscMatrix &Q, Rcpp::List prec) {
  Q.p = Rcpp::as<std::vector<int>>(prec["Q_p"]);
  Q.i = Rcpp::as<std::vector<int>>(prec["Q_i"]);
  if (Q.p.empty()) {
    Rcpp::stop("random_fem Q_p must be non-empty");
  }
  Q.nrow_ = static_cast<int>(Q.p.size()) - 1;
  Q.ncol_ = Q.nrow_;
  Q.x.clear();
}

void fill_fem_L1_pattern(CscMatrix &L1, Rcpp::List chol) {
  if (!chol.containsElementNamed("perm") || !chol.containsElementNamed("L1")) {
    Rcpp::stop("random_fem chol must contain perm and L1");
  }
  L1 = CscMatrix(Rcpp::S4(chol["L1"]));
}

FemCholPayload read_fem_payload(const density_data &model) {
  if (Rf_isNull(model.precision) || TYPEOF(model.precision) != VECSXP) {
    Rcpp::stop(
        "random_fem precision must be a list from fem_precision_payload()");
  }
  Rcpp::List prec(model.precision);
  for (const char *key :
       {"Q_p", "Q_i", "C_x", "G_x", "G2_x", "chol", "alpha"}) {
    if (!prec.containsElementNamed(key)) {
      Rcpp::stop("random_fem precision missing '%s'", key);
    }
  }
  if (model.num_theta < 2) {
    Rcpp::stop("random_fem requires two theta parameters (range, sd)");
  }
  FemCholPayload out;
  out.alpha = Rcpp::as<int>(prec["alpha"]);
  fill_fem_Q_pattern(out.Q, prec);
  out.C_x = Rcpp::as<std::vector<double>>(prec["C_x"]);
  out.G_x = Rcpp::as<std::vector<double>>(prec["G_x"]);
  out.G2_x = Rcpp::as<std::vector<double>>(prec["G2_x"]);
  out.G3_x = prec.containsElementNamed("G3_x")
                 ? Rcpp::as<std::vector<double>>(prec["G3_x"])
                 : std::vector<double>(out.Q.nnz(), 0.0);
  Rcpp::List chol = prec["chol"];
  out.perm = Rcpp::as<std::vector<int>>(chol["perm"]);
  fill_fem_L1_pattern(out.L1, chol);
  if (out.C_x.size() != out.Q.nnz() || out.G_x.size() != out.Q.nnz() ||
      out.G2_x.size() != out.Q.nnz() || out.G3_x.size() != out.Q.nnz()) {
    Rcpp::stop("random_fem Gram coefficient vectors must match length(Q_i)");
  }
  return out;
}

// Symmetry weights for the upper-triangle storage: 1 on the diagonal, 2 off
// it, so a single pass over the stored nonzeros gives a full quadratic form.
std::vector<double> fem_symmetry_weights(const FemCholPayload &pay) {
  std::vector<double> w(pay.Q.nnz());
  for (std::size_t col = 0; col + 1 < pay.Q.p.size(); ++col) {
    for (int pos = pay.Q.p[col]; pos < pay.Q.p[col + 1]; ++pos) {
      const std::size_t row =
          static_cast<std::size_t>(pay.Q.i[static_cast<std::size_t>(pos)]);
      w[static_cast<std::size_t>(pos)] = (row == col) ? 1.0 : 2.0;
    }
  }
  return w;
}

// Register the constant part of Q(theta) with the fem_logdet atomic:
// Q = a*C + b*G + c*G2 (+ d*G3), so the atomic only sees the coefficient
// scalars while the Grams, pattern, and symbolic factor stay fixed.
template <int Alpha>
std::size_t register_fem_logdet_payload(const FemCholPayload &pay) {
  femlogdet::Payload atom;
  atom.Q = pay.Q;
  atom.n = static_cast<std::size_t>(pay.Q.ncol());
  atom.w = fem_symmetry_weights(pay);
  atom.M.push_back(pay.C_x);
  atom.M.push_back(pay.G_x);
  atom.M.push_back(pay.G2_x);
  if (Alpha == 3) {
    atom.M.push_back(pay.G3_x);
  }
  atom.perm = pay.perm;
  atom.perm_inv.assign(atom.n, 0);
  for (std::size_t r = 0; r < atom.n; ++r) {
    atom.perm_inv[static_cast<std::size_t>(pay.perm[r])] = static_cast<int>(r);
  }
  atom.L1 = pay.L1;
  return femlogdet::fem_logdet_atomic_instance().register_payload(
      std::move(atom));
}

// Which tape slots hold the two thetas, and whether each already arrives
// logged (transform_theta) or has to be logged here.
struct FemThetaSpec {
  std::size_t t0_tape = 0;
  std::size_t t1_tape = 0;
  bool range_is_log = false;
  bool sd_is_log = false;
};

FemThetaSpec fem_theta_spec(const density_data &model, const Config &config) {
  FemThetaSpec spec;
  spec.t0_tape = model.theta_index(0);
  spec.t1_tape = model.theta_index(1);
  spec.range_is_log = transform_theta_at(config, model.theta_row(0));
  spec.sd_is_log = transform_theta_at(config, model.theta_row(1));
  return spec;
}

// Coefficients c of Q(theta) = sum_j c_j * M_j with M = (C, G, G2, [G3]).
//
// Public theta is (range, sd): practical range rho = sqrt(8*nu)/kappa and
// field SD, with nu = Alpha - 1 in 2D. This is the single definition of the
// theta -> coefficient map; the quadratic form, the log-determinant atomic,
// and the analytic shard's coefficient tape all go through it so they cannot
// drift apart.
template <int Alpha>
void fem_Q_coefficients(const CppAD::AD<double> &raw_range,
                        const CppAD::AD<double> &raw_sd, bool range_is_log,
                        bool sd_is_log,
                        CppAD::vector<CppAD::AD<double>> &out) {
  CppAD::AD<double> log_range = raw_range;
  CppAD::AD<double> log_sd = raw_sd;
  if (!range_is_log) {
    log_range = CppAD::log(log_range);
  }
  if (!sd_is_log) {
    log_sd = CppAD::log(log_sd);
  }
  const CppAD::AD<double> range = CppAD::exp(log_range);
  const CppAD::AD<double> sd = CppAD::exp(log_sd);
  const CppAD::AD<double> kappa =
      CppAD::sqrt(CppAD::AD<double>(8.0 * (Alpha - 1))) / range;
  const CppAD::AD<double> k2 = kappa * kappa;
  const CppAD::AD<double> k4 = k2 * k2;
  // Alpha=2: sigma^2 = 1/(4 pi kappa^2 tau^2)
  // Alpha=3: sigma^2 = 1/(8 pi kappa^4 tau^2)
  const CppAD::AD<double> tau =
      (Alpha == 2) ? CppAD::AD<double>(1) /
                         (kappa * sd * CppAD::sqrt(CppAD::AD<double>(PIx4)))
                   : CppAD::AD<double>(1) /
                         (k2 * sd * CppAD::sqrt(CppAD::AD<double>(PIx8)));
  const CppAD::AD<double> tau2 = tau * tau;

  out.resize(Alpha == 2 ? 3 : 4);
  if constexpr (Alpha == 2) {
    out[0] = tau2 * k4;
    out[1] = CppAD::AD<double>(2) * tau2 * k2;
    out[2] = tau2;
  } else {
    const CppAD::AD<double> k6 = k4 * k2;
    out[0] = tau2 * k6;
    out[1] = CppAD::AD<double>(3) * tau2 * k4;
    out[2] = CppAD::AD<double>(3) * tau2 * k2;
    out[3] = tau2;
  }
}

// Domain 2, Range m tape of the coefficient map, used by the analytic shard
// to get c, dc/dtheta, and d2c/dtheta2 without hand-differentiating the
// range/sd -> kappa/tau composition or the optional log transform.
template <int Alpha>
CppAD::ADFun<double> make_fem_coef_fun(const FemThetaSpec &spec, double seed0,
                                       double seed1) {
  if (!spec.range_is_log && !(seed0 > 0.0)) {
    seed0 = 1.0;
  }
  if (!spec.sd_is_log && !(seed1 > 0.0)) {
    seed1 = 1.0;
  }
  CppAD::vector<CppAD::AD<double>> ax(2);
  ax[0] = seed0;
  ax[1] = seed1;
  CppAD::Independent(ax);
  CppAD::vector<CppAD::AD<double>> ay;
  fem_Q_coefficients<Alpha>(ax[0], ax[1], spec.range_is_log, spec.sd_is_log,
                            ay);
  return CppAD::ADFun<double>(ax, ay);
}

// Owned CSC for one Gram: shared Q pattern + that Gram's value vector.
CscMatrix gram_csc(const CscMatrix &pattern, const std::vector<double> &vals) {
  CscMatrix m;
  m.p = pattern.p;
  m.i = pattern.i;
  m.x = vals;
  m.nrow_ = pattern.nrow_;
  m.ncol_ = pattern.ncol_;
  return m;
}

// Taped path: qf = sum_j c_j(theta) * (gamma' M_j gamma) via quadform atomics
// so the tape records O(m) atomic calls instead of O(nnz(Q)) products.
template <int Alpha>
CppAD::vector<CppAD::AD<double>>
random_fem_ssq(const CppAD::vector<CppAD::AD<double>> &x, const density_data &model,
               const Config &config) {
  adlaplace::quadform::init_quadform_atomic();
  const FemCholPayload pay = read_fem_payload(model);
  if (pay.alpha != Alpha) {
    Rcpp::stop("precision alpha (%d) does not match kernel alpha (%d)",
               pay.alpha, Alpha);
  }
  const FemThetaSpec spec = fem_theta_spec(model, config);
  CppAD::vector<CppAD::AD<double>> c;
  fem_Q_coefficients<Alpha>(x[spec.t0_tape], x[spec.t1_tape],
                            spec.range_is_log, spec.sd_is_log, c);

  const std::size_t n = static_cast<std::size_t>(pay.Q.ncol());
  const std::vector<std::size_t> gidx = model.all_gamma_global_indices();
  if (gidx.size() != n) {
    Rcpp::stop("random_fem: length(gamma) (%d) != nrow(Q) (%d)",
               static_cast<int>(gidx.size()), static_cast<int>(n));
  }

  CppAD::vector<CppAD::AD<double>> ax(n);
  for (std::size_t j = 0; j < n; ++j) {
    ax[j] = x[gidx[j]];
  }

  const std::size_t m = c.size();
  const std::vector<double> *grams[4] = {&pay.C_x, &pay.G_x, &pay.G2_x,
                                         &pay.G3_x};
  CppAD::AD<double> qf = 0.0;
  for (std::size_t j = 0; j < m; ++j) {
    const std::size_t call_id = adlaplace::quadform::register_csc(
        gram_csc(pay.Q, *grams[j]));
    CppAD::vector<CppAD::AD<double>> ay(1);
    adlaplace::quadform::call_quadform(call_id, ax, ay);
    qf += c[j] * ay[0];
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = -CppAD::AD<double>(0.5) * qf;
  return result;
}

template <int Alpha>
CppAD::vector<CppAD::AD<double>>
random_fem_det(const CppAD::vector<CppAD::AD<double>> &x, const density_data &model,
               const Config &config) {
  const FemCholPayload pay = read_fem_payload(model);
  if (pay.alpha != Alpha) {
    Rcpp::stop("precision alpha (%d) does not match kernel alpha (%d)",
               pay.alpha, Alpha);
  }
  const std::size_t n = static_cast<std::size_t>(pay.Q.ncol());
  const std::size_t call_id = register_fem_logdet_payload<Alpha>(pay);

  // Coefficients of Q = a*C + b*G + c*G2 (+ d*G3); the atomic handles the
  // factorization and Tr(Q^{-1} dQ) derivatives in doubles, off the tape.
  const FemThetaSpec spec = fem_theta_spec(model, config);
  CppAD::vector<CppAD::AD<double>> ax;
  fem_Q_coefficients<Alpha>(x[spec.t0_tape], x[spec.t1_tape],
                            spec.range_is_log, spec.sd_is_log, ax);
  CppAD::vector<CppAD::AD<double>> ay(1);
  femlogdet::fem_logdet_atomic_instance()(call_id, ax, ay);
  const CppAD::AD<double> log_det = ay[0];

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] =
      CppAD::AD<double>(0.5) * log_det - CppAD::AD<double>(n * ONEHALFLOGTWOPI);
  return result;
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_ssq_sparsity(const density_data &model) {
  const FemCholPayload pay = read_fem_payload(model);
  const std::size_t idx_range = model.theta_index(0);
  const std::size_t idx_sd = model.theta_index(1);
  const std::vector<std::size_t> gidx = model.all_gamma_global_indices();

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  for (std::size_t col = 0; col < gidx.size(); ++col) {
    for (int pos = pay.Q.p[col]; pos < pay.Q.p[col + 1]; ++pos) {
      const std::size_t row =
          static_cast<std::size_t>(pay.Q.i[static_cast<std::size_t>(pos)]);
      pairs.insert({gidx[row], gidx[col]});
      pairs.insert({gidx[col], gidx[row]});
    }
  }
  for (std::size_t k = 0; k < gidx.size(); ++k) {
    pairs.insert({gidx[k], idx_range});
    pairs.insert({idx_range, gidx[k]});
    pairs.insert({gidx[k], idx_sd});
    pairs.insert({idx_sd, gidx[k]});
  }
  pairs.insert({idx_range, idx_range});
  pairs.insert({idx_sd, idx_sd});
  pairs.insert({idx_range, idx_sd});
  pairs.insert({idx_sd, idx_range});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, model.n_tape, pairs);
  return hessian;
}

// Gradient sparsity of -0.5 gamma' Q(theta) gamma: every gamma tape slot
// (Q has a full diagonal) plus the two thetas. Columns are emitted in
// ascending order to match for_jac_sparsity.
CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
fem_ssq_grad_sparsity(const density_data &model) {
  std::set<std::size_t> cols;
  for (std::size_t g : model.all_gamma_global_indices()) {
    cols.insert(g);
  }
  cols.insert(model.theta_index(0));
  cols.insert(model.theta_index(1));

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad;
  grad.resize(1, model.n_tape, cols.size());
  std::size_t k = 0;
  for (std::size_t j : cols) {
    grad.set(k++, 0, j);
  }
  return grad;
}

template <int Alpha>
AdTape build_fem_ssq_pack(density_data model, const Config &cfg,
                          bool verbose) {
  validate_config_matches_model(cfg, model, false);
  if (Rf_isNull(model.precision)) {
    Rcpp::stop("precision is required for random densities");
  }
  model.apply_tape_domain(cfg, "all", 0);
  AdTape pack;
  pack.owner_thread_assigned = false;
  adpack_attach_tape_maps(pack, model);
  AdpackPatterns pats = adpack_build_patterns(
      pack, model.n_tape, model.seq_gamma, fem_ssq_grad_sparsity(model),
      random_fem_ssq_sparsity(model), verbose);
  adpack_install_patterns(pack, pats);
  return pack;
}

CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
random_fem_det_sparsity(const density_data &model) {
  const std::size_t idx_range = model.theta_index(0);
  const std::size_t idx_sd = model.theta_index(1);

  std::set<std::pair<std::size_t, std::size_t>> pairs;
  pairs.insert({idx_range, idx_range});
  pairs.insert({idx_sd, idx_sd});
  pairs.insert({idx_range, idx_sd});
  pairs.insert({idx_sd, idx_range});

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  set_sparse_rc_pairs(hessian, model.n_tape, pairs);
  return hessian;
}

// Immutable data for the analytic ssq shard. Built from the compacted model
// so gidx / theta indices are tape positions, matching the tape the patterns
// were discovered on.
template <int Alpha>
std::shared_ptr<const femssq::Data>
make_fem_ssq_data(density_data model, const Config &cfg) {
  model.apply_tape_domain(cfg, "all", 0);
  const FemCholPayload pay = read_fem_payload(model);

  auto data = std::make_shared<femssq::Data>();
  data->alpha = Alpha;
  data->n = static_cast<std::size_t>(pay.Q.ncol());
  data->Q = pay.Q;
  data->w = fem_symmetry_weights(pay);
  data->M.push_back(pay.C_x);
  data->M.push_back(pay.G_x);
  data->M.push_back(pay.G2_x);
  if (Alpha == 3) {
    data->M.push_back(pay.G3_x);
  }

  data->gidx = model.all_gamma_global_indices();
  if (data->gidx.size() != data->n) {
    Rcpp::stop("random_fem: length(gamma) (%d) != nrow(Q) (%d)",
               static_cast<int>(data->gidx.size()),
               static_cast<int>(data->n));
  }

  const FemThetaSpec spec = fem_theta_spec(model, cfg);
  data->t_tape[0] = spec.t0_tape;
  data->t_tape[1] = spec.t1_tape;
  const double seed0 = cfg.theta.size() > 0 ? cfg.theta[0] : 1.0;
  const double seed1 = cfg.theta.size() > 1 ? cfg.theta[1] : 1.0;
  data->coef_fun = make_fem_coef_fun<Alpha>(spec, seed0, seed1);
  return data;
}

template <int Alpha>
SEXP create_ad_shard_random_fem_ssq(SEXP model, Rcpp::List config) {
  const density_data ad_model(model);
  const Config cfg(config);
  // Two modes only:
  //   fem_analytic=TRUE  (default)  untaped FemSsqShard closed forms
  //   fem_analytic=FALSE            EvalShard over a quadform-atomic tape
  const bool fem_analytic = adlaplace_get_bool(config, "fem_analytic", true);

  std::vector<AdTape> packs;
  if (fem_analytic) {
    packs.push_back(build_fem_ssq_pack<Alpha>(ad_model, cfg, cfg.verbose));
    return make_ad_pack_ptr(fem_ssq_packs_to_ad_fun(
        std::move(packs), ad_model.num_beta, ad_model.num_theta,
        make_fem_ssq_data<Alpha>(ad_model, cfg)));
  }
  packs.push_back(build_ad_fun_random_with_pattern(
      ad_model, config, random_fem_ssq<Alpha>, random_fem_ssq_sparsity));
  return make_ad_pack_ptr(packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                          ad_model.num_theta,
                                          adlaplace_fem_make_shard));
}

LogDensSingleDataFn resolve_fem_det(const std::string &name) {
  if (name == "random_fem_det_2") {
    return random_fem_det<2>;
  }
  if (name == "random_fem_det_3") {
    return random_fem_det<3>;
  }
  Rcpp::stop("unknown parameters density: %s", name.c_str());
}

// Numeric LDL log-determinant of Q(range, sd) for debugging. Returns the
// same path the random_fem_det_* atomic uses, plus the D pivots.
template <int Alpha>
Rcpp::List fem_logdet_debug_impl(const FemCholPayload &pay, double range,
                                 double sd, bool range_is_log,
                                 bool sd_is_log) {
  if (pay.alpha != Alpha) {
    Rcpp::stop("precision alpha (%d) does not match requested alpha (%d)",
               pay.alpha, Alpha);
  }
  FemThetaSpec spec;
  spec.t0_tape = 0;
  spec.t1_tape = 1;
  spec.range_is_log = range_is_log;
  spec.sd_is_log = sd_is_log;

  CppAD::vector<CppAD::AD<double>> ax(2);
  ax[0] = range;
  ax[1] = sd;
  CppAD::Independent(ax);
  CppAD::vector<CppAD::AD<double>> ay;
  fem_Q_coefficients<Alpha>(ax[0], ax[1], spec.range_is_log, spec.sd_is_log,
                            ay);
  CppAD::ADFun<double> coef_fun(ax, ay);

  CPPAD_TESTVECTOR(double) tv(2);
  tv[0] = range;
  tv[1] = sd;
  const CPPAD_TESTVECTOR(double) cv = coef_fun.Forward(0, tv);
  const std::size_t m = cv.size();
  Rcpp::NumericVector coef(m);
  for (std::size_t j = 0; j < m; ++j) {
    coef[j] = cv[j];
  }

  const std::size_t nnz = pay.Q.nnz();
  const std::size_t n = static_cast<std::size_t>(pay.Q.ncol());
  std::vector<double> Q_x(nnz, 0.0);
  for (std::size_t k = 0; k < nnz; ++k) {
    double v = coef[0] * pay.C_x[k] + coef[1] * pay.G_x[k] +
               coef[2] * pay.G2_x[k];
    if constexpr (Alpha == 3) {
      v += coef[3] * pay.G3_x[k];
    }
    Q_x[k] = v;
  }

  std::vector<double> L_x(pay.L1.nnz(), 0.0);
  std::vector<double> D(n, 0.0);
  const double log_det = adlaplace::chol::chol_update_csc(
      pay.Q.p, pay.Q.i, Q_x, pay.perm, pay.L1.p, pay.L1.i, L_x, D);

  int n_bad = 0;
  int first_bad = NA_INTEGER;
  for (std::size_t j = 0; j < n; ++j) {
    if (D[j] <= 0.0 || !std::isfinite(D[j])) {
      if (n_bad == 0) {
        first_bad = static_cast<int>(j);
      }
      ++n_bad;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("logdet") = log_det, Rcpp::Named("coef") = coef,
      Rcpp::Named("D") = Rcpp::NumericVector(D.begin(), D.end()),
      Rcpp::Named("first_bad") = first_bad, Rcpp::Named("n_bad") = n_bad,
      Rcpp::Named("n") = static_cast<int>(n), Rcpp::Named("alpha") = Alpha);
}

FemCholPayload read_fem_payload_list(Rcpp::List prec) {
  for (const char *key :
       {"Q_p", "Q_i", "C_x", "G_x", "G2_x", "chol", "alpha"}) {
    if (!prec.containsElementNamed(key)) {
      Rcpp::stop("fem_logdet_debug: precision missing '%s'", key);
    }
  }
  FemCholPayload out;
  out.alpha = Rcpp::as<int>(prec["alpha"]);
  fill_fem_Q_pattern(out.Q, prec);
  out.C_x = Rcpp::as<std::vector<double>>(prec["C_x"]);
  out.G_x = Rcpp::as<std::vector<double>>(prec["G_x"]);
  out.G2_x = Rcpp::as<std::vector<double>>(prec["G2_x"]);
  out.G3_x = prec.containsElementNamed("G3_x")
                 ? Rcpp::as<std::vector<double>>(prec["G3_x"])
                 : std::vector<double>(out.Q.nnz(), 0.0);
  Rcpp::List chol = prec["chol"];
  out.perm = Rcpp::as<std::vector<int>>(chol["perm"]);
  fill_fem_L1_pattern(out.L1, chol);
  return out;
}

} // namespace

//' Sparse LDL log-determinant of a FEM precision at given (range, sd)
//'
//' Exercises the same \code{chol_update_csc} path used by
//' \code{random_fem_det_*}. Useful for diagnosing non-finite joint densities.
//'
//' @param precision List from \code{\link{fem_precision_payload}}.
//' @param range,sd Practical range and field SD (natural or log scale).
//' @param log_scale Logical length-1 or length-2: whether \code{range}/\code{sd}
//'   are already on the log scale. Default \code{FALSE} (\code{NULL}).
//' @return List with \code{logdet}, \code{coef}, \code{D}, \code{first_bad},
//'   \code{n_bad}, \code{n}, \code{alpha}.
//' @export
// [[Rcpp::export]]
Rcpp::List fem_logdet_debug(
    Rcpp::List precision, double range, double sd,
    Rcpp::Nullable<Rcpp::LogicalVector> log_scale = R_NilValue) {
  const FemCholPayload pay = read_fem_payload_list(precision);
  bool range_is_log = false;
  bool sd_is_log = false;
  if (log_scale.isNotNull()) {
    Rcpp::LogicalVector ls(log_scale);
    if (ls.size() >= 1) {
      range_is_log = ls[0];
    }
    if (ls.size() >= 2) {
      sd_is_log = ls[1];
    } else if (ls.size() == 1) {
      sd_is_log = ls[0];
    }
  }
  if (pay.alpha == 2) {
    return fem_logdet_debug_impl<2>(pay, range, sd, range_is_log, sd_is_log);
  }
  if (pay.alpha == 3) {
    return fem_logdet_debug_impl<3>(pay, range, sd, range_is_log, sd_is_log);
  }
  Rcpp::stop("fem_logdet_debug: alpha must be 2 or 3, got %d", pay.alpha);
}

//' Build raw AD handle for a random_fem_ssq_2 term
//'
//' @param model An \code{density_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP create_ad_shard_random_fem_ssq_2(SEXP model, Rcpp::List config) {
  return create_ad_shard_random_fem_ssq<2>(model, config);
}

//' Build raw AD handle for a random_fem_ssq_3 term
//'
//' @param model An \code{density_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP create_ad_shard_random_fem_ssq_3(SEXP model, Rcpp::List config) {
  return create_ad_shard_random_fem_ssq<3>(model, config);
}

//' Build parameters-shard \code{ad_pack_ptr} for FEM log-determinant densities.
//'
//' @param model \code{density_data} S4 object with FEM precision payload.
//' @param config Model configuration list.
//' @param name Density name (\code{"random_fem_det_2"} or \code{"random_fem_det_3"}).
//' @return External pointer of class \code{ad_pack_ptr}.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
SEXP get_ad_pack_raw_parameters(SEXP model, Rcpp::List config, std::string name) {
  density_data ad_model(model);
  const Config cfg(config);
  ad_model.apply_tape_domain(cfg, "all", 0);
  AdTape pack = build_ad_fun_parameters(
      ad_model, config, resolve_fem_det(name),
      random_fem_det_sparsity(ad_model));
  std::vector<AdTape> packs;
  packs.push_back(std::move(pack));
  return make_ad_pack_ptr(packs_to_ad_fun(std::move(packs), ad_model.num_beta,
                                         ad_model.num_theta,
                                         adlaplace_fem_make_shard));
}
