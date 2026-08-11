#include "sun32_tape.hpp"

#include "mvn_sparsity.hpp"
#include "ompad.hpp"
#include "pmvn_atomic.hpp"
#include "qnorm_atomic.hpp"
#include "sun_ad_utils.hpp"

#include <Rcpp.h>

#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace admvn {
namespace {

constexpr double kLog2Pi = 1.8378770664093453;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

using AD = CppAD::AD<double>;
using ChAD = std::vector<std::vector<AD>>;
using Mat = std::vector<std::vector<double>>;

ChAD zeros_ad(std::size_t n) {
  return ChAD(n, std::vector<AD>(n, AD(0.0)));
}

ChAD zeros_ad_dm() {
  return ChAD(kSun32D, std::vector<AD>(kSun32M, AD(0.0)));
}

void chol_lower_ad(const ChAD& a, ChAD& l) {
  const std::size_t n = a.size();
  l = zeros_ad(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      AD s = a[i][j];
      for (std::size_t k = 0; k < j; ++k) {
        s -= l[i][k] * l[j][k];
      }
      if (i == j) {
        l[i][j] = CppAD::sqrt(CppAD::CondExpGt(s, AD(1e-12), s, AD(1e-12)));
      } else {
        l[i][j] = s / l[j][j];
      }
    }
  }
}

ChAD cov2cor(const ChAD& sigma) {
  const std::size_t n = sigma.size();
  std::vector<AD> scale(n);
  for (std::size_t i = 0; i < n; ++i) {
    const AD s2 = CppAD::CondExpGt(sigma[i][i], AD(1e-12), sigma[i][i], AD(1e-12));
    scale[i] = CppAD::sqrt(s2);
  }
  ChAD out = zeros_ad(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      out[i][j] = sigma[i][j] / (scale[i] * scale[j]);
    }
  }
  return out;
}

// Unit-row Cholesky for joint dim = D + M = 5; 10 free coordinates.
ChAD unit_row_chol_joint(const std::vector<AD>& z) {
  if (z.size() != kSun32HsNFree) {
    throw std::runtime_error("unit_row_chol_joint: expected 10 free parameters");
  }
  ChAD Br = zeros_ad(kSun32JointD);
  Br[0][0] = AD(1.0);
  std::size_t idx = 0;
  const AD kEps = AD(1e-6);
  for (std::size_t i = 1; i < kSun32JointD; ++i) {
    AD sum_sq = AD(0.0);
    for (std::size_t j = 0; j < i; ++j) {
      const AD tij = CppAD::tanh(z[idx++]);
      const AD rem = CppAD::CondExpGt(
        AD(1.0) - sum_sq, kEps, AD(1.0) - sum_sq, kEps);
      Br[i][j] = tij * CppAD::sqrt(rem);
      sum_sq += Br[i][j] * Br[i][j];
    }
    const AD rem_diag = CppAD::CondExpGt(
      AD(1.0) - sum_sq, kEps, AD(1.0) - sum_sq, kEps);
    Br[i][i] = CppAD::sqrt(rem_diag);
  }
  return Br;
}

struct Sun32DpAD {
  std::vector<AD> xi = std::vector<AD>(kSun32D, AD(0.0));
  ChAD Omega = zeros_ad(kSun32D);          // D x D
  ChAD Delta = zeros_ad_dm();              // D x M
  ChAD Gamma = zeros_ad(kSun32M);          // M x M
};

// Parameter order (length 16): xi(3), nu(3), z(10)
Sun32DpAD make_sun32_dp_hs_ad(const std::vector<AD>& par) {
  if (par.size() != kSun32NPar) {
    throw std::runtime_error("par must have length 16");
  }
  Sun32DpAD out;
  for (std::size_t i = 0; i < kSun32D; ++i) {
    out.xi[i] = par[i];
  }
  std::vector<AD> z(kSun32HsNFree);
  for (std::size_t i = 0; i < kSun32HsNFree; ++i) {
    z[i] = par[2 * kSun32D + i];
  }
  ChAD B = unit_row_chol_joint(z);
  ChAD J = zeros_ad(kSun32JointD);
  for (std::size_t i = 0; i < kSun32JointD; ++i) {
    for (std::size_t j = 0; j < kSun32JointD; ++j) {
      for (std::size_t k = 0; k < kSun32JointD; ++k) {
        J[i][j] += B[i][k] * B[j][k];
      }
    }
  }

  // Omega = D_u C_u D_u from U-block of J
  for (std::size_t i = 0; i < kSun32D; ++i) {
    const AD si = CppAD::CondExpGt(
      par[kSun32D + i], AD(1e-12), par[kSun32D + i], AD(1e-12));
    for (std::size_t j = 0; j < kSun32D; ++j) {
      const AD sj = CppAD::CondExpGt(
        par[kSun32D + j], AD(1e-12), par[kSun32D + j], AD(1e-12));
      out.Omega[i][j] = si * sj * J[i][j];
    }
  }
  // Delta = J_UV (D x M)
  for (std::size_t i = 0; i < kSun32D; ++i) {
    for (std::size_t j = 0; j < kSun32M; ++j) {
      out.Delta[i][j] = J[i][kSun32D + j];
    }
  }
  // Gamma = J_VV (M x M)
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32M; ++j) {
      out.Gamma[i][j] = J[kSun32D + i][kSun32D + j];
    }
  }
  return out;
}

Sun32DpAD make_sun32_dp_dispatch(
  const std::vector<AD>& par, Sun32ParMap par_map) {
  if (par_map == Sun32ParMap::kHyperspherical) {
    return make_sun32_dp_hs_ad(par);
  }
  throw std::runtime_error(
    "SUN(3,2) block-Cholesky parameterization is not implemented; use hs");
}

Mat zeros_mat(std::size_t n) {
  return Mat(n, std::vector<double>(n, 0.0));
}

// Returns Omega (D x D); fills lambda (M x M).
Mat sun32_sigma_double(
  const std::vector<double>& par, Mat& lambda, Sun32ParMap par_map) {
  std::vector<AD> par_ad(kSun32NPar);
  for (std::size_t i = 0; i < kSun32NPar; ++i) {
    par_ad[i] = AD(par[i]);
  }
  Sun32DpAD dp = make_sun32_dp_dispatch(par_ad, par_map);

  ChAD omega_bar = cov2cor(dp.Omega);
  ChAD omega_inv = zeros_ad(kSun32D);
  ChAD l{};
  chol_lower_ad(omega_bar, l);
  for (std::size_t col = 0; col < kSun32D; ++col) {
    std::vector<AD> y(kSun32D, AD(0.0));
    for (std::size_t i = 0; i < kSun32D; ++i) {
      y[i] = (i == col) ? AD(1.0) : AD(0.0);
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      AD s = y[i];
      for (std::size_t k = 0; k < i; ++k) {
        s -= l[i][k] * y[k];
      }
      y[i] = s / l[i][i];
    }
    for (int i = static_cast<int>(kSun32D) - 1; i >= 0; --i) {
      AD s = y[static_cast<std::size_t>(i)];
      for (std::size_t k = static_cast<std::size_t>(i) + 1; k < kSun32D; ++k) {
        s -= l[k][static_cast<std::size_t>(i)] * y[k];
      }
      y[static_cast<std::size_t>(i)] =
        s / l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      omega_inv[i][col] = y[i];
    }
  }

  // d_oinv = Delta' Omega^{-1}  (M x D)
  ChAD d_oinv(kSun32M, std::vector<AD>(kSun32D, AD(0.0)));
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32D; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSun32D; ++k) {
        s += dp.Delta[k][i] * omega_inv[k][j];
      }
      d_oinv[i][j] = s;
    }
  }

  // lambda = Gamma - Delta' Omega^{-1} Delta  (M x M)
  ChAD lambda_ad = dp.Gamma;
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32M; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSun32D; ++k) {
        s += d_oinv[i][k] * dp.Delta[k][j];
      }
      lambda_ad[i][j] -= s;
    }
  }

  lambda = zeros_mat(kSun32M);
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32M; ++j) {
      lambda[i][j] = Value(lambda_ad[i][j]);
    }
  }
  Mat omega = zeros_mat(kSun32D);
  for (std::size_t i = 0; i < kSun32D; ++i) {
    for (std::size_t j = 0; j < kSun32D; ++j) {
      omega[i][j] = Value(dp.Omega[i][j]);
    }
  }
  return omega;
}

AD log_mvn_obs(const std::vector<AD>& x, const Sun32DpAD& dp) {
  ChAD l{};
  chol_lower_ad(dp.Omega, l);
  AD log_det = AD(0.0);
  for (std::size_t i = 0; i < kSun32D; ++i) {
    log_det += CppAD::log(l[i][i]);
  }
  log_det *= AD(2.0);

  std::vector<AD> y(kSun32D);
  for (std::size_t i = 0; i < kSun32D; ++i) {
    y[i] = x[i] - dp.xi[i];
  }
  for (std::size_t i = 0; i < kSun32D; ++i) {
    AD s = y[i];
    for (std::size_t k = 0; k < i; ++k) {
      s -= l[i][k] * y[k];
    }
    y[i] = s / l[i][i];
  }
  AD quad = AD(0.0);
  for (std::size_t i = 0; i < kSun32D; ++i) {
    quad += y[i] * y[i];
  }
  return AD(-0.5 * static_cast<double>(kSun32D)) * AD(kLog2Pi)
    - AD(0.5) * log_det - AD(0.5) * quad;
}

ChAD lambda_from_dp(const Sun32DpAD& dp) {
  ChAD omega_bar = cov2cor(dp.Omega);
  ChAD omega_inv = zeros_ad(kSun32D);
  ChAD l{};
  chol_lower_ad(omega_bar, l);
  for (std::size_t col = 0; col < kSun32D; ++col) {
    std::vector<AD> y(kSun32D, AD(0.0));
    for (std::size_t i = 0; i < kSun32D; ++i) {
      y[i] = (i == col) ? AD(1.0) : AD(0.0);
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      AD s = y[i];
      for (std::size_t k = 0; k < i; ++k) {
        s -= l[i][k] * y[k];
      }
      y[i] = s / l[i][i];
    }
    for (int i = static_cast<int>(kSun32D) - 1; i >= 0; --i) {
      AD s = y[static_cast<std::size_t>(i)];
      for (std::size_t k = static_cast<std::size_t>(i) + 1; k < kSun32D; ++k) {
        s -= l[k][static_cast<std::size_t>(i)] * y[k];
      }
      y[static_cast<std::size_t>(i)] =
        s / l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      omega_inv[i][col] = y[i];
    }
  }

  ChAD d_oinv(kSun32M, std::vector<AD>(kSun32D, AD(0.0)));
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32D; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSun32D; ++k) {
        s += dp.Delta[k][i] * omega_inv[k][j];
      }
      d_oinv[i][j] = s;
    }
  }

  ChAD lambda = dp.Gamma;
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32M; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSun32D; ++k) {
        s += d_oinv[i][k] * dp.Delta[k][j];
      }
      lambda[i][j] -= s;
    }
  }
  return lambda;
}

// alpha has length M (latent orthant dimension)
std::vector<AD> alpha_from_x(const std::vector<AD>& x, const Sun32DpAD& dp) {
  ChAD omega_bar = cov2cor(dp.Omega);
  ChAD omega_inv = zeros_ad(kSun32D);
  ChAD l{};
  chol_lower_ad(omega_bar, l);
  for (std::size_t col = 0; col < kSun32D; ++col) {
    std::vector<AD> y(kSun32D, AD(0.0));
    for (std::size_t i = 0; i < kSun32D; ++i) {
      y[i] = (i == col) ? AD(1.0) : AD(0.0);
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      AD s = y[i];
      for (std::size_t k = 0; k < i; ++k) {
        s -= l[i][k] * y[k];
      }
      y[i] = s / l[i][i];
    }
    for (int i = static_cast<int>(kSun32D) - 1; i >= 0; --i) {
      AD s = y[static_cast<std::size_t>(i)];
      for (std::size_t k = static_cast<std::size_t>(i) + 1; k < kSun32D; ++k) {
        s -= l[k][static_cast<std::size_t>(i)] * y[k];
      }
      y[static_cast<std::size_t>(i)] =
        s / l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      omega_inv[i][col] = y[i];
    }
  }

  std::vector<AD> z(kSun32D);
  for (std::size_t i = 0; i < kSun32D; ++i) {
    const AD omega_i = CppAD::sqrt(dp.Omega[i][i]);
    z[i] = (x[i] - dp.xi[i]) / omega_i;
  }

  ChAD d_oinv(kSun32M, std::vector<AD>(kSun32D, AD(0.0)));
  for (std::size_t i = 0; i < kSun32M; ++i) {
    for (std::size_t j = 0; j < kSun32D; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSun32D; ++k) {
        s += dp.Delta[k][i] * omega_inv[k][j];
      }
      d_oinv[i][j] = s;
    }
  }

  std::vector<AD> alpha(kSun32M);
  for (std::size_t i = 0; i < kSun32M; ++i) {
    AD s = AD(0.0);
    for (std::size_t j = 0; j < kSun32D; ++j) {
      s += d_oinv[i][j] * z[j];
    }
    alpha[i] = s;
  }
  return alpha;
}

// Pack GenZ inputs for an n-dimensional covariance (n = M for lambda/Gamma).
void pack_genz_ad(
  const ChAD& sigma,
  std::vector<AD>& scale,
  std::vector<AD>& vech_ch) {

  const std::size_t n = sigma.size();
  ChAD corr = cov2cor(sigma);
  ChAD ch{};
  chol_lower_ad(corr, ch);
  scale.assign(n, AD(0.0));
  for (std::size_t i = 0; i < n; ++i) {
    const AD s2 = CppAD::CondExpGt(sigma[i][i], AD(1e-12), sigma[i][i], AD(1e-12));
    scale[i] = CppAD::sqrt(s2);
  }
  vech_ch.assign(vech_size(n), AD(0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      vech_ch[vech_index(i, j)] = ch[i][j];
    }
  }
}

AD sun32_obs_contrib_ad(
  const std::vector<AD>& x,
  const Sun32DpAD& dp,
  const ChAD& /*lambda*/,
  const std::vector<AD>& scale_l,
  const std::vector<AD>& vech_l,
  MvnTape& p1_tape,
  std::size_t pmvn_call_id) {

  const AD log_pdf = log_mvn_obs(x, dp);
  const std::vector<AD> alpha = alpha_from_x(x, dp);
  std::vector<AD> mean_v(kSun32M, AD(0.0));
  const AD p1 = detail::pmvn_ad(
    pmvn_call_id, p1_tape, alpha, mean_v, scale_l, vech_l);
  return log_pdf + log_prob_floored(p1);
}

AD sun32_obs_loglik_ad(
  const std::vector<AD>& x,
  const std::vector<AD>& par,
  MvnTape& p1_tape,
  std::size_t pmvn_call_id,
  Sun32ParMap par_map) {

  Sun32DpAD dp = make_sun32_dp_dispatch(par, par_map);
  const ChAD lambda = lambda_from_dp(dp);
  std::vector<AD> scale_l;
  std::vector<AD> vech_l;
  pack_genz_ad(lambda, scale_l, vech_l);
  return sun32_obs_contrib_ad(
    x, dp, lambda, scale_l, vech_l, p1_tape, pmvn_call_id);
}

AD sun32_log_p2_ad(
  const std::vector<AD>& par,
  MvnTape& p2_tape,
  std::size_t pmvn_call_id,
  Sun32ParMap par_map) {

  Sun32DpAD dp = make_sun32_dp_dispatch(par, par_map);
  std::vector<AD> scale_g;
  std::vector<AD> vech_g;
  pack_genz_ad(dp.Gamma, scale_g, vech_g);

  std::vector<AD> mean_v(kSun32M, AD(0.0));
  std::vector<AD> upper_zero(kSun32M, AD(0.0));
  const AD p2 = detail::pmvn_ad(
    pmvn_call_id, p2_tape, upper_zero, mean_v, scale_g, vech_g);
  return log_prob_floored(p2);
}

void setup_adfun_sparsity(Sun32AdfunPack& pack, const std::vector<double>& par_seed) {
  pack.x_seed = detail::to_cppad_vector(par_seed);
  std::vector<std::size_t> par_subset(kSun32NPar);
  std::iota(par_subset.begin(), par_subset.end(), 0);
  detail::setup_mvn_sparsity(
    pack.fun,
    pack.x_seed,
    par_subset,
    pack.pattern_grad,
    pack.pattern_grad_inner,
    pack.pattern_hessian,
    pack.pattern_hessian_inner,
    pack.work_grad,
    pack.work_inner_grad,
    pack.work_hess,
    pack.work_inner_hess,
    pack.unused_pattern,
    pack.w);
}

void build_obs_shard_fun(
  Sun32ObsShard& shard,
  const std::vector<double>& par_seed,
  Sun32ParMap par_map) {
  detail::init_qnorm_atomic();
  detail::init_pmvn_atomic();
  detail::set_pmvn_atomic_tape(shard.pmvn_call_id, &shard.p1_tape);

  std::vector<AD> x_ad(kSun32D);
  for (std::size_t i = 0; i < kSun32D; ++i) {
    x_ad[i] = shard.x[i];
  }

  std::vector<AD> par_ad(kSun32NPar);
  for (std::size_t i = 0; i < kSun32NPar; ++i) {
    par_ad[i] = AD(par_seed[i]);
  }
  CppAD::Independent(par_ad);

  AD y = sun32_obs_loglik_ad(
    x_ad, par_ad, shard.p1_tape, shard.pmvn_call_id, par_map);
  std::vector<AD> out(1);
  out[0] = y;
  shard.adfun.fun = CppAD::ADFun<double>(par_ad, out);
}

void build_p2_shard_fun(
  Sun32P2Shard& shard,
  const std::vector<double>& par_seed,
  Sun32ParMap par_map) {
  detail::init_qnorm_atomic();
  detail::init_pmvn_atomic();
  detail::set_pmvn_atomic_tape(shard.pmvn_call_id, &shard.p2_tape);

  std::vector<AD> par_ad(kSun32NPar);
  for (std::size_t i = 0; i < kSun32NPar; ++i) {
    par_ad[i] = AD(par_seed[i]);
  }
  CppAD::Independent(par_ad);

  AD y = sun32_log_p2_ad(
    par_ad, shard.p2_tape, shard.pmvn_call_id, par_map);
  std::vector<AD> out(1);
  out[0] = y;
  shard.adfun.fun = CppAD::ADFun<double>(par_ad, out);
}

Sun32Result eval_adfun_log(
  Sun32AdfunPack& pack,
  std::size_t pmvn_call_id,
  MvnTape* pmvn_tape,
  const std::vector<double>& par,
  int deriv) {

  Sun32Result out;
  out.gradient.assign(kSun32NPar, 0.0);
  out.hessian.assign(kSun32NPar, std::vector<double>(kSun32NPar, 0.0));

  detail::init_qnorm_atomic();
  detail::init_pmvn_atomic();
  detail::set_pmvn_atomic_tape(pmvn_call_id, pmvn_tape);

  const CppAD::vector<double> u = detail::to_cppad_vector(par);
  const double raw = pack.fun.Forward(0, u)[0];
  out.value = raw;
  out.error = 0.0;

  if (deriv < 1) {
    return out;
  }

  CppAD::vector<double> w_rev(1);
  w_rev[0] = 1.0;
  const CppAD::vector<double> dw = pack.fun.Reverse(1, w_rev);
  for (std::size_t j = 0; j < kSun32NPar; ++j) {
    out.gradient[j] += dw[j];
  }

  if (deriv < 2) {
    return out;
  }

  CppAD::sparse_rc<CppAD::vector<size_t>> hes_pattern;
  pack.fun.sparse_hes(u, pack.w, pack.pattern_hessian, hes_pattern, detail::kHessColor, pack.work_hess);

  const auto& rows = pack.pattern_hessian.row();
  const auto& hcols = pack.pattern_hessian.col();
  const auto& hvals = pack.pattern_hessian.val();
  for (size_t k = 0; k < pack.pattern_hessian.nnz(); ++k) {
    const std::size_t i = rows[k];
    const std::size_t j = hcols[k];
    out.hessian[i][j] += hvals[k];
    if (i != j) {
      out.hessian[j][i] += hvals[k];
    }
  }

  return out;
}

void apply_log_scale(Sun32Result& out, bool log_scale, int deriv) {
  if (log_scale) {
    return;
  }
  const std::vector<double> grad_log = out.gradient;
  const std::vector<std::vector<double>> h_log = out.hessian;
  const double f = std::exp(out.value);
  out.value = f;
  if (deriv >= 1) {
    for (std::size_t i = 0; i < kSun32NPar; ++i) {
      out.gradient[i] = f * grad_log[i];
    }
  }
  if (deriv >= 2) {
    for (std::size_t i = 0; i < kSun32NPar; ++i) {
      for (std::size_t j = 0; j < kSun32NPar; ++j) {
        out.hessian[i][j] = f * (h_log[i][j] + grad_log[i] * grad_log[j]);
      }
    }
  }
}

void accumulate_result(
  Sun32Result& total,
  const Sun32Result& part,
  double scale,
  int deriv) {

  total.value += scale * part.value;
  if (deriv >= 1) {
    for (std::size_t i = 0; i < kSun32NPar; ++i) {
      total.gradient[i] += scale * part.gradient[i];
    }
  }
  if (deriv >= 2) {
    for (std::size_t i = 0; i < kSun32NPar; ++i) {
      for (std::size_t j = 0; j < kSun32NPar; ++j) {
        total.hessian[i][j] += scale * part.hessian[i][j];
      }
    }
  }
}

int resolve_n_threads(int n_threads, int bundle_default) {
  if (n_threads <= 0) {
    n_threads = bundle_default;
  }
  if (n_threads < 1) {
    n_threads = 1;
  }
#ifdef _OPENMP
  const int max_procs = omp_get_num_procs();
  if (max_procs > 0 && n_threads > max_procs) {
    n_threads = max_procs;
  }
#else
  n_threads = 1;
#endif
  return n_threads;
}

}  // namespace

Sun32TapeBundle create_sun32_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights,
  Sun32ParMap par_map) {

  if (x_rows.empty() || par_seed.size() != kSun32NPar) {
    throw std::runtime_error("invalid seed dimensions for SUN(3,2) tape");
  }
  for (const auto& row : x_rows) {
    if (row.size() != kSun32D) {
      throw std::runtime_error("each observation must have length 3");
    }
  }
  if (!weights.empty() && weights.size() != x_rows.size()) {
    throw std::runtime_error("weights must have length equal to number of rows in x");
  }

  Sun32TapeBundle bundle;
  bundle.n_obs = x_rows.size();
  bundle.n_threads = resolve_n_threads(n_threads, 1);
  bundle.par_map = par_map;
  if (weights.empty()) {
    bundle.weights.assign(x_rows.size(), 1.0);
  } else {
    bundle.weights = weights;
  }
  bundle.weight_sum = 0.0;
  for (double w : bundle.weights) {
    bundle.weight_sum += w;
  }

  Mat lambda = zeros_mat(kSun32M);
  (void)sun32_sigma_double(par_seed, lambda, par_map);
  Mat gamma = zeros_mat(kSun32M);
  {
    std::vector<AD> par_ad(kSun32NPar);
    for (std::size_t i = 0; i < kSun32NPar; ++i) {
      par_ad[i] = AD(par_seed[i]);
    }
    Sun32DpAD dp = make_sun32_dp_dispatch(par_ad, par_map);
    for (std::size_t i = 0; i < kSun32M; ++i) {
      for (std::size_t j = 0; j < kSun32M; ++j) {
        gamma[i][j] = Value(dp.Gamma[i][j]);
      }
    }
  }

  // Latent orthant is M-dimensional
  std::vector<double> lower(kSun32M, kNegInf);
  std::vector<double> mean_zero(kSun32M, 0.0);

  const std::size_t p2_call_id = x_rows.size();
  bundle.p2.pmvn_call_id = p2_call_id;

  bundle.p2.p2_tape = create_mvn_tape(
    lower, mean_zero, mean_zero, gamma,
    n_points, n_shifts, seed + 1U, /*value_only=*/false, /*reorder=*/false);
  build_p2_shard_fun(bundle.p2, par_seed, par_map);
  setup_adfun_sparsity(bundle.p2.adfun, par_seed);

  bundle.shards.reserve(x_rows.size());
  for (std::size_t obs = 0; obs < x_rows.size(); ++obs) {
    Sun32ObsShard shard;
    shard.x = x_rows[obs];
    shard.pmvn_call_id = obs;

    std::vector<AD> par_ad(kSun32NPar);
    std::vector<AD> x_ad(kSun32D);
    for (std::size_t i = 0; i < kSun32NPar; ++i) {
      par_ad[i] = AD(par_seed[i]);
    }
    for (std::size_t i = 0; i < kSun32D; ++i) {
      x_ad[i] = AD(x_rows[obs][i]);
    }
    Sun32DpAD dp = make_sun32_dp_dispatch(par_ad, par_map);
    std::vector<AD> alpha_ad = alpha_from_x(x_ad, dp);
    std::vector<double> upper_alpha(kSun32M);
    for (std::size_t i = 0; i < kSun32M; ++i) {
      upper_alpha[i] = Value(alpha_ad[i]);
    }

    shard.p1_tape = create_mvn_tape(
      lower, upper_alpha, mean_zero, lambda,
      n_points, n_shifts, seed + 2U + static_cast<unsigned int>(obs),
      /*value_only=*/false, /*reorder=*/false);

    build_obs_shard_fun(shard, par_seed, par_map);
    setup_adfun_sparsity(shard.adfun, par_seed);
    bundle.shards.push_back(std::move(shard));
  }

  return bundle;
}

Sun32Result eval_sun32_bundle(
  Sun32TapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale,
  int deriv,
  int n_threads,
  bool manage_parallel_scope) {

  if (par.size() != kSun32NPar) {
    throw std::runtime_error("par must have length 16");
  }
  if (deriv < 0 || deriv > 2) {
    throw std::runtime_error("deriv must be 0, 1, or 2");
  }

  n_threads = resolve_n_threads(n_threads, bundle.n_threads);

  Sun32Result total;
  total.gradient.assign(kSun32NPar, 0.0);
  total.hessian.assign(kSun32NPar, std::vector<double>(kSun32NPar, 0.0));

  std::unique_ptr<CppadParallelScope> parallel_scope;
  if (manage_parallel_scope) {
    parallel_scope =
      std::make_unique<CppadParallelScope>(static_cast<std::size_t>(n_threads));
  }

  const std::size_t n_shards = bundle.shards.size();

#ifdef _OPENMP
  if (n_threads > 1 && n_shards > 0) {
    std::vector<double> shard_values(n_shards, 0.0);
    std::vector<std::vector<double>> shard_grads;
    std::vector<std::vector<std::vector<double>>> shard_hess;
    if (deriv >= 1) {
      shard_grads.assign(n_shards, std::vector<double>(kSun32NPar, 0.0));
    }
    if (deriv >= 2) {
      shard_hess.assign(
        n_shards,
        std::vector<std::vector<double>>(kSun32NPar, std::vector<double>(kSun32NPar, 0.0)));
    }

#pragma omp parallel num_threads(n_threads)
    {
#pragma omp for schedule(static)
      for (std::size_t i = 0; i < n_shards; ++i) {
        Sun32Result shard_res = eval_adfun_log(
          bundle.shards[i].adfun,
          bundle.shards[i].pmvn_call_id,
          &bundle.shards[i].p1_tape,
          par,
          deriv);
        shard_values[i] = shard_res.value;
        if (deriv >= 1) {
          shard_grads[i] = std::move(shard_res.gradient);
        }
        if (deriv >= 2) {
          shard_hess[i] = std::move(shard_res.hessian);
        }
      }
    }

    for (std::size_t i = 0; i < n_shards; ++i) {
      const double wi = bundle.weights[i];
      total.value += wi * shard_values[i];
      if (deriv >= 1) {
        for (std::size_t j = 0; j < kSun32NPar; ++j) {
          total.gradient[j] += wi * shard_grads[i][j];
        }
      }
      if (deriv >= 2) {
        for (std::size_t r = 0; r < kSun32NPar; ++r) {
          for (std::size_t c = 0; c < kSun32NPar; ++c) {
            total.hessian[r][c] += wi * shard_hess[i][r][c];
          }
        }
      }
    }
  } else
#endif
  {
    for (std::size_t i = 0; i < n_shards; ++i) {
      Sun32Result shard_res = eval_adfun_log(
        bundle.shards[i].adfun,
        bundle.shards[i].pmvn_call_id,
        &bundle.shards[i].p1_tape,
        par,
        deriv);
      accumulate_result(total, shard_res, bundle.weights[i], deriv);
    }
  }

  Sun32Result p2_res = eval_adfun_log(
    bundle.p2.adfun,
    bundle.p2.pmvn_call_id,
    &bundle.p2.p2_tape,
    par,
    deriv);
  accumulate_result(total, p2_res, -bundle.weight_sum, deriv);

  apply_log_scale(total, log_scale, deriv);
  return total;
}

}  // namespace admvn
