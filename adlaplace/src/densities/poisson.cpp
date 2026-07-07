#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/densities/poisson.hpp"

#include <cmath>

// Poisson observation density with an optional offset (e.g. log expected
// counts for standardized-mortality models):
//   y_i ~ Poisson(exp(eta_i + offset_i)),  eta = X beta + A gamma.
//
// The offset is read from config$offset (numeric, length ny or scalar;
// missing/NULL means zero). There is no theta parameter. Terms that do not
// depend on (beta, gamma) -- y_i * offset_i and -lgamma(y_i + 1) -- are
// accumulated in plain double and enter the tape as constants, so no
// separate "extra" parameters shard is needed.
CppAD::vector<CppAD::AD<double>> poisson_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup) {

  const Config config(config_list);

  Rcpp::NumericVector offset;
  bool have_offset = config_list.containsElementNamed("offset") &&
    !Rf_isNull(config_list["offset"]);
  if (have_offset) {
    offset = Rcpp::as<Rcpp::NumericVector>(config_list["offset"]);
    if (offset.size() != 1 &&
        static_cast<size_t>(offset.size()) != model.y.size()) {
      Rcpp::stop(
        "config$offset must have length 1 or length(y) (%d), got %d",
        static_cast<int>(model.y.size()), static_cast<int>(offset.size()));
    }
  }

  const bool have_shards = config.shards.ncol() > 0;
  const size_t ny = model.y.size();
  size_t startP = 0;
  size_t endP = 0;
  if (have_shards) {
    startP = config.shards.p[Dgroup];
    endP = config.shards.p[Dgroup + 1];
  } else if (Dgroup == 0) {
    endP = ny;
  }

  CppAD::AD<double> linear = 0.0;
  CppAD::AD<double> expsum = 0.0;
  double constants = 0.0;

  for (size_t DI = startP; DI < endP; ++DI) {
    const size_t Dobs = have_shards ? config.shards.i[DI] : DI;

    CppAD::AD<double> eta = 0.0;
    const size_t p0x = model.XTp.p[Dobs];
    const size_t p1x = model.XTp.p[Dobs + 1];
    for (size_t D = p0x; D < p1x; ++D) {
      eta += model.XTp.x[D] * x[model.XTp.i[D]];
    }
    const size_t p0a = model.ATp.p[Dobs];
    const size_t p1a = model.ATp.p[Dobs + 1];
    for (size_t D = p0a; D < p1a; ++D) {
      eta += model.ATp.x[D] * x[model.num_beta + model.ATp.i[D]];
    }

    const double off = have_offset
      ? (offset.size() == 1 ? offset[0] : offset[Dobs])
      : 0.0;
    const double y = model.y[Dobs];

    linear += y * eta;
    expsum += CppAD::exp(eta + off);
    constants += y * off - std::lgamma(y + 1.0);
  }

  CppAD::vector<CppAD::AD<double>> result(1);
  result[0] = linear - expsum + CppAD::AD<double>(constants);
  return result;
}
