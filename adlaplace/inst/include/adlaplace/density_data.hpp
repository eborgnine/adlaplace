#ifndef ADLAPLACE_AD_DATA_HPP
#define ADLAPLACE_AD_DATA_HPP

#include "adlaplace/rviews.hpp"

#include <Rinternals.h>
#include <vector>
#include <cstddef>

struct ad_data {
  DgCView beta_map;
  DgCView gamma_map;
  DgCView theta_map;

  DgCView ATp;
  DgCView XTp;
  DgCView elgm_matrix;
  NumVecView y;
  NumVecView weights;

  SEXP precision = R_NilValue;

  std::size_t num_beta = 0;
  std::size_t num_gamma = 0;
  std::size_t num_theta = 0;
  std::size_t num_full = 0;

  // Global x indices: local design column j -> x[beta_global[j]], etc.
  std::vector<std::size_t> beta_global;
  std::vector<std::size_t> gamma_global;
  std::vector<std::size_t> theta_global;

  // Global parameter indices for gamma (beta + 0..num_gamma-1).
  std::vector<int> seq_gamma;

  explicit ad_data(SEXP data_sexp);

  std::size_t theta_index(int col = 0) const;
  std::vector<std::size_t> gamma_global_indices(int col = 0) const;
  std::vector<std::size_t> all_gamma_global_indices() const;
  DgCView mult_precision_Q() const;
  double mult_precision_rank() const;
  double mult_precision_log_det() const;
};

inline SEXP ad_data_slot(SEXP obj, const char* name) {
  return R_do_slot(obj, Rf_mkString(name));
}

namespace adlaplace_detail {

inline std::vector<std::size_t> build_map_global_lookup(
  const DgCView& map,
  std::size_t block_offset,
  std::size_t empty_sentinel,
  const char* what,
  bool require_nonempty_columns) {

  const int ncols = map.ncol();
  std::vector<std::size_t> out(static_cast<std::size_t>(ncols), empty_sentinel);
  for (int col = 0; col < ncols; ++col) {
    const int start = map.p[col];
    const int end = map.p[col + 1];
    const int nnz_col = end - start;
    if (nnz_col == 0) {
      if (require_nonempty_columns) {
        Rcpp::stop(
          "%s column %d must have exactly one nonzero (got none)",
          what,
          col + 1
        );
      }
      continue;
    }
    if (nnz_col != 1) {
      Rcpp::stop(
        "%s column %d must have exactly one nonzero (got %d)",
        what,
        col + 1,
        nnz_col
      );
    }
    const int local_row = map.i[start];
    if (local_row < 0 || local_row >= map.nrow()) {
      Rcpp::stop(
        "%s column %d has row index %d out of range [0, %d)",
        what,
        col + 1,
        local_row,
        map.nrow()
      );
    }
    out[static_cast<std::size_t>(col)] =
      block_offset + static_cast<std::size_t>(local_row);
  }
  return out;
}

}  // namespace adlaplace_detail

inline ad_data::ad_data(SEXP data_sexp) {
  if (!Rf_inherits(data_sexp, "ad_data")) {
    Rcpp::stop("expected S4 object of class ad_data");
  }

  beta_map = DgCView(Rcpp::as<Rcpp::S4>(ad_data_slot(data_sexp, "beta_map")));
  gamma_map = DgCView(Rcpp::as<Rcpp::S4>(ad_data_slot(data_sexp, "gamma_map")));
  theta_map = DgCView(Rcpp::as<Rcpp::S4>(ad_data_slot(data_sexp, "theta_map")));

  ATp = DgCView(Rcpp::as<Rcpp::S4>(ad_data_slot(data_sexp, "ATp")));
  XTp = DgCView(Rcpp::as<Rcpp::S4>(ad_data_slot(data_sexp, "XTp")));
  elgm_matrix = DgCView(Rcpp::as<Rcpp::S4>(ad_data_slot(data_sexp, "elgm_matrix")));
  y = NumVecView(ad_data_slot(data_sexp, "y"));
  weights = NumVecView(ad_data_slot(data_sexp, "weights"));
  precision = ad_data_slot(data_sexp, "precision");

  num_beta = static_cast<std::size_t>(beta_map.nrow());
  num_gamma = static_cast<std::size_t>(gamma_map.nrow());
  num_theta = static_cast<std::size_t>(theta_map.nrow());
  num_full = num_beta + num_gamma + num_theta;

  const std::size_t theta_sentinel = num_beta + num_gamma;
  const bool have_beta_design = XTp.nrow() > 0;
  const bool have_gamma_design = ATp.nrow() > 0;

  if (have_beta_design && beta_map.ncol() != XTp.nrow()) {
    Rcpp::stop(
      "beta_map ncol (%d) must match ncol(X) (%d)",
      beta_map.ncol(),
      XTp.nrow()
    );
  }
  if (have_gamma_design && gamma_map.ncol() != ATp.nrow()) {
    Rcpp::stop(
      "gamma_map ncol (%d) must match ncol(A) (%d)",
      gamma_map.ncol(),
      ATp.nrow()
    );
  }

  beta_global = adlaplace_detail::build_map_global_lookup(
    beta_map,
    0,
    num_full,
    "beta_map",
    have_beta_design
  );
  gamma_global = adlaplace_detail::build_map_global_lookup(
    gamma_map,
    num_beta,
    num_full,
    "gamma_map",
    have_gamma_design
  );
  theta_global = adlaplace_detail::build_map_global_lookup(
    theta_map,
    num_beta + num_gamma,
    theta_sentinel,
    "theta_map",
    false
  );

  seq_gamma.resize(num_gamma);
  for (std::size_t d = 0; d < num_gamma; ++d) {
    seq_gamma[d] = static_cast<int>(num_beta + d);
  }
}

inline std::size_t ad_data::theta_index(int col) const {
  if (col < 0 || col >= static_cast<int>(theta_global.size())) {
    return num_beta + num_gamma;
  }
  return theta_global[static_cast<std::size_t>(col)];
}

inline std::vector<std::size_t> ad_data::gamma_global_indices(int col) const {
  std::vector<std::size_t> out;
  if (col < 0 || col >= static_cast<int>(gamma_global.size())) {
    return out;
  }
  const std::size_t g = gamma_global[static_cast<std::size_t>(col)];
  if (g >= num_full) {
    return out;
  }
  out.push_back(g);
  return out;
}

inline std::vector<std::size_t> ad_data::all_gamma_global_indices() const {
  std::vector<std::size_t> out;
  out.reserve(gamma_global.size());
  for (std::size_t g : gamma_global) {
    if (g < num_full) {
      out.push_back(g);
    }
  }
  return out;
}

inline DgCView ad_data::mult_precision_Q() const {
  if (Rf_isNull(precision) || TYPEOF(precision) != VECSXP) {
    Rcpp::stop(
        "random_mult precision must be list(Q = <dgCMatrix>, log_det, rank)");
  }
  Rcpp::List prec(precision);
  if (!prec.containsElementNamed("Q")) {
    Rcpp::stop("random_mult precision list must contain Q");
  }
  SEXP q_sexp = prec["Q"];
  if (!Rf_inherits(q_sexp, "dgCMatrix") && !Rf_inherits(q_sexp, "ngCMatrix")) {
    Rcpp::stop("random_mult precision Q must be a dgCMatrix "
               "(convert with as(Q, \"generalMatrix\") for symmetric storage)");
  }
  return DgCView(Rcpp::as<Rcpp::S4>(q_sexp));
}

inline double ad_data::mult_precision_rank() const {
  const int n_term = gamma_map.ncol();
  if (Rf_isNull(precision) || TYPEOF(precision) != VECSXP) {
    return static_cast<double>(n_term);
  }
  Rcpp::List prec(precision);
  return prec.containsElementNamed("rank")
    ? Rcpp::as<double>(prec["rank"])
    : static_cast<double>(n_term);
}

inline double ad_data::mult_precision_log_det() const {
  if (Rf_isNull(precision) || TYPEOF(precision) != VECSXP) {
    return 0.0;
  }
  Rcpp::List prec(precision);
  return prec.containsElementNamed("log_det")
    ? Rcpp::as<double>(prec["log_det"])
    : 0.0;
}

inline void validate_config_matches_model(
  const Config& cfg,
  const ad_data& data,
  bool check_gamma = true) {
  if (cfg.beta.size() != data.num_beta) {
    Rcpp::stop(
      "length(config$beta) (%d) must match nrow(beta_map) (%d)",
      static_cast<int>(cfg.beta.size()),
      static_cast<int>(data.num_beta)
    );
  }
  if (check_gamma && cfg.gamma.size() != data.num_gamma) {
    Rcpp::stop(
      "length(config$gamma) (%d) must match nrow(gamma_map) (%d)",
      static_cast<int>(cfg.gamma.size()),
      static_cast<int>(data.num_gamma)
    );
  }
  if (cfg.theta.size() != data.num_theta) {
    Rcpp::stop(
      "length(config$theta) (%d) must match nrow(theta_map) (%d)",
      static_cast<int>(cfg.theta.size()),
      static_cast<int>(data.num_theta)
    );
  }
}

#endif
