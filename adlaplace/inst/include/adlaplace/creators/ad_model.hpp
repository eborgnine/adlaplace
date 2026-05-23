#ifndef ADLAPLACE_AD_MODEL_HPP
#define ADLAPLACE_AD_MODEL_HPP

#include "adlaplace/creators/rviews.hpp"

#include <Rinternals.h>
#include <vector>
#include <cstddef>

struct ad_model {
  DgCView beta_map;
  DgCView gamma_map;
  DgCView theta_map;

  DgCView ATp;
  DgCView XTp;
  DgCView elgm_matrix;
  NumVecView y;

  std::size_t num_beta = 0;
  std::size_t num_gamma = 0;
  std::size_t num_theta = 0;
  std::size_t num_parameters = 0;
  std::size_t num_full = 0;

  // Global parameter indices for gamma (beta + 0..num_gamma-1).
  std::vector<int> seq_gamma;

  explicit ad_model(SEXP model_sexp);

  std::size_t theta_index(int col = 0) const;
  std::vector<std::size_t> gamma_global_indices(int col = 0) const;
  std::vector<std::size_t> all_gamma_global_indices() const;
};

inline SEXP ad_model_slot(SEXP obj, const char* name) {
  return R_do_slot(obj, Rf_mkString(name));
}

inline ad_model::ad_model(SEXP model_sexp) {
  if (!Rf_inherits(model_sexp, "ad_model")) {
    Rcpp::stop("expected S4 object of class ad_model");
  }

  beta_map = DgCView(Rcpp::as<Rcpp::S4>(ad_model_slot(model_sexp, "beta_map")));
  gamma_map = DgCView(Rcpp::as<Rcpp::S4>(ad_model_slot(model_sexp, "gamma_map")));
  theta_map = DgCView(Rcpp::as<Rcpp::S4>(ad_model_slot(model_sexp, "theta_map")));

  ATp = DgCView(Rcpp::as<Rcpp::S4>(ad_model_slot(model_sexp, "ATp")));
  XTp = DgCView(Rcpp::as<Rcpp::S4>(ad_model_slot(model_sexp, "XTp")));
  elgm_matrix = DgCView(Rcpp::as<Rcpp::S4>(ad_model_slot(model_sexp, "elgmMatrix")));
  y = NumVecView(ad_model_slot(model_sexp, "y"));

  num_beta = static_cast<std::size_t>(beta_map.nrow());
  num_gamma = static_cast<std::size_t>(gamma_map.nrow());
  num_theta = static_cast<std::size_t>(theta_map.nrow());
  num_parameters = num_beta + num_theta;
  num_full = num_beta + num_gamma + num_theta;

  seq_gamma.resize(num_gamma);
  for (std::size_t d = 0; d < num_gamma; ++d) {
    seq_gamma[d] = static_cast<int>(num_beta + d);
  }
}

inline std::size_t ad_model::theta_index(int col) const {
  if (col < 0 || col >= theta_map.ncol()) {
    return num_beta + num_gamma;
  }
  const int start = theta_map.p[col];
  const int end = theta_map.p[col + 1];
  if (start >= end) {
    return num_beta + num_gamma;
  }
  const int local_row = theta_map.i[start];
  return num_beta + num_gamma + static_cast<std::size_t>(local_row);
}

inline std::vector<std::size_t> ad_model::gamma_global_indices(int col) const {
  std::vector<std::size_t> out;
  if (col < 0 || col >= gamma_map.ncol()) {
    return out;
  }
  const int start = gamma_map.p[col];
  const int end = gamma_map.p[col + 1];
  out.reserve(static_cast<std::size_t>(end - start));
  for (int k = start; k < end; ++k) {
    out.push_back(num_beta + static_cast<std::size_t>(gamma_map.i[k]));
  }
  return out;
}

inline std::vector<std::size_t> ad_model::all_gamma_global_indices() const {
  std::vector<std::size_t> out;
  out.reserve(static_cast<std::size_t>(gamma_map.ncol()));
  for (int col = 0; col < gamma_map.ncol(); ++col) {
    const std::vector<std::size_t> idx = gamma_global_indices(col);
    out.insert(out.end(), idx.begin(), idx.end());
  }
  return out;
}

inline void validate_config_matches_model(
  const Config& cfg,
  const ad_model& model,
  bool check_gamma = true) {
  if (cfg.beta.size() != model.num_beta) {
    Rcpp::stop(
      "length(config$beta) (%d) must match nrow(beta_map) (%d)",
      static_cast<int>(cfg.beta.size()),
      static_cast<int>(model.num_beta)
    );
  }
  if (check_gamma && cfg.gamma.size() != model.num_gamma) {
    Rcpp::stop(
      "length(config$gamma) (%d) must match nrow(gamma_map) (%d)",
      static_cast<int>(cfg.gamma.size()),
      static_cast<int>(model.num_gamma)
    );
  }
  if (cfg.theta.size() != model.num_theta) {
    Rcpp::stop(
      "length(config$theta) (%d) must match nrow(theta_map) (%d)",
      static_cast<int>(cfg.theta.size()),
      static_cast<int>(model.num_theta)
    );
  }
}

#endif
