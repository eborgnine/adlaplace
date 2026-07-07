#include "adlaplace/creators/ad_data.hpp"
#include "adlaplace/densities/binomial.hpp"

// Bernoulli (binomial, size 1) observation density with a logit link:
//   y_i ~ Bernoulli(logit^{-1}(eta_i)),  eta = X beta + A gamma.
// Per-point log density: y_i * eta_i - log(1 + exp(eta_i)).
// There is no theta parameter and no separate extra shard.
CppAD::vector<CppAD::AD<double>> binomial_obs(
  const CppAD::vector<CppAD::AD<double>>& x,
  const ad_data& model,
  const Rcpp::List& config_list,
  const size_t Dgroup) {

  const Config config(config_list);

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

  CppAD::AD<double> result = 0.0;
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

    const CppAD::AD<double> softplus = CppAD::CondExpGt(
      eta, CppAD::AD<double>(0.0),
      eta + CppAD::log1p(CppAD::exp(-eta)),
      CppAD::log1p(CppAD::exp(eta)));

    result += CppAD::AD<double>(model.y[Dobs]) * eta - softplus;
  }

  CppAD::vector<CppAD::AD<double>> out(1);
  out[0] = result;
  return out;
}
