#include "sun_tape.hpp"

#include "mvn_sparsity.hpp"
#include "ompad.hpp"
#include "pmvn_atomic.hpp"
#include "qnorm_atomic.hpp"

#include <Rcpp.h>

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace admvn {
namespace {

constexpr double kLog2Pi = 1.8378770664093453;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

using AD = CppAD::AD<double>;
using ChAD = std::vector<std::vector<AD>>;
using Mat3 = std::array<std::array<double, kSunD>, kSunD>;

ChAD zeros_ad3() {
  return ChAD(kSunD, std::vector<AD>(kSunD, AD(0.0)));
}

ChAD identity_ad3() {
  ChAD m = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    m[i][i] = AD(1.0);
  }
  return m;
}

void chol_lower_ad3(const ChAD& a, ChAD& l) {
  l = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
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

ChAD mat_mul3(const ChAD& a, const ChAD& b) {
  ChAD out = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        s += a[i][k] * b[k][j];
      }
      out[i][j] = s;
    }
  }
  return out;
}

ChAD mat_transpose3(const ChAD& a) {
  ChAD out = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      out[i][j] = a[j][i];
    }
  }
  return out;
}

ChAD diag3(const std::array<AD, kSunD>& d) {
  ChAD out = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    out[i][i] = d[i];
  }
  return out;
}

ChAD cov2cor3(const ChAD& sigma) {
  std::array<AD, kSunD> scale{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    const AD s2 = CppAD::CondExpGt(sigma[i][i], AD(1e-12), sigma[i][i], AD(1e-12));
    scale[i] = CppAD::sqrt(s2);
  }
  ChAD out = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      out[i][j] = sigma[i][j] / (scale[i] * scale[j]);
    }
  }
  return out;
}

// Unit-row Cholesky (hyperspherical) factor of a correlation matrix.
  // Free params z (length d*(d-1)/2) are unconstrained; each partial coordinate
  // uses tanh so C = Br Br^T is a residual correlation (unit diagonal).
ChAD unit_row_chol3(const std::array<AD, 3>& z) {
  ChAD Br = zeros_ad3();
  Br[0][0] = AD(1.0);
  std::size_t idx = 0;
  for (std::size_t i = 1; i < kSunD; ++i) {
    AD sum_sq = AD(0.0);
    for (std::size_t j = 0; j < i; ++j) {
      const AD tij = CppAD::tanh(z[idx++]);
      const AD rem = CppAD::CondExpGt(AD(1.0) - sum_sq, AD(1e-16),
                                      AD(1.0) - sum_sq, AD(1e-16));
      Br[i][j] = tij * CppAD::sqrt(rem);
      sum_sq += Br[i][j] * Br[i][j];
    }
    const AD rem_diag = CppAD::CondExpGt(AD(1.0) - sum_sq, AD(1e-16),
                                         AD(1.0) - sum_sq, AD(1e-16));
    Br[i][i] = CppAD::sqrt(rem_diag);
  }
  return Br;
}

struct SunDpAD {
  std::array<AD, kSunD> xi{};
  ChAD Omega = zeros_ad3();
  ChAD Delta = zeros_ad3();
  ChAD Gamma = zeros_ad3();
};

SunDpAD make_sun_dp_ad(const std::vector<AD>& par) {
  SunDpAD out;
  for (std::size_t i = 0; i < kSunD; ++i) {
    out.xi[i] = par[i];
  }

  // Correlation C_u from ell; scale S = diag(nu) applied outside.
  // L acts on C_u so |L_ij| <= 1 is the natural regime.
  const std::array<AD, kSunD> nu{par[3], par[4], par[5]};
  const std::array<AD, 3> z_ell{par[6], par[7], par[8]};
  ChAD Bu = unit_row_chol3(z_ell);
  ChAD Cu = mat_mul3(Bu, mat_transpose3(Bu));
  std::array<AD, kSunD> s_u{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    // nu = marginal standard deviations
    s_u[i] = CppAD::CondExpGt(nu[i], AD(1e-12), nu[i], AD(1e-12));
  }
  ChAD Omega = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      Omega[i][j] = s_u[i] * s_u[j] * Cu[i][j];
    }
  }
  out.Omega = Omega;

  // L: diagonals; pair levels L12,L13,L23; signed gaps e12,e13,e23
  // lower = Lij + eij, upper = Lij - eij
  const AD L11 = par[9];
  const AD L22 = par[10];
  const AD L33 = par[11];
  const AD L12 = par[12];
  const AD L13 = par[13];
  const AD L23 = par[14];
  const AD e12 = par[15];
  const AD e13 = par[16];
  const AD e23 = par[17];
  ChAD L = zeros_ad3();
  L[0][0] = L11;
  L[0][1] = L12 - e12;
  L[0][2] = L13 - e13;
  L[1][0] = L12 + e12;
  L[1][1] = L22;
  L[1][2] = L23 - e23;
  L[2][0] = L13 + e13;
  L[2][1] = L23 + e23;
  L[2][2] = L33;

  const std::array<AD, 3> z_br{par[18], par[19], par[20]};
  ChAD Br = unit_row_chol3(z_br);
  ChAD C = mat_mul3(Br, mat_transpose3(Br));

  // Correlation-scale block: M = L C_u L'; Delta = C_u L'
  ChAD M = mat_mul3(mat_mul3(L, Cu), mat_transpose3(L));

  // Residual R = D^{1/2} C D^{1/2} with D = diag(1 - diag(M)), so
  // Sigma_vv = M + R has unit diagonal (correlation matrix).
  std::array<AD, kSunD> s_res{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    const AD rem = CppAD::CondExpGt(AD(1.0) - M[i][i], AD(1e-12),
                                    AD(1.0) - M[i][i], AD(1e-12));
    s_res[i] = CppAD::sqrt(rem);
  }
  ChAD R = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      R[i][j] = s_res[i] * s_res[j] * C[i][j];
    }
  }

  ChAD Sigma_vv = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      Sigma_vv[i][j] = M[i][j] + R[i][j];
    }
  }

  // Delta = D_u^{-1} Sigma_uv = C_u L'  (S cancels)
  out.Delta = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD acc = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        acc += Cu[i][k] * L[j][k];
      }
      out.Delta[i][j] = acc;
    }
  }

  out.Gamma = Sigma_vv;
  return out;
}

Mat3 sun_sigma_double(const std::vector<double>& par, Mat3& lambda) {
  std::vector<AD> par_ad(kSunNPar);
  for (std::size_t i = 0; i < kSunNPar; ++i) {
    par_ad[i] = AD(par[i]);
  }
  SunDpAD dp = make_sun_dp_ad(par_ad);

  ChAD omega_bar = cov2cor3(dp.Omega);
  ChAD omega_inv = zeros_ad3();
  ChAD l{};
  chol_lower_ad3(omega_bar, l);
  for (std::size_t col = 0; col < kSunD; ++col) {
    std::array<AD, kSunD> y{};
    for (std::size_t i = 0; i < kSunD; ++i) {
      y[i] = (i == col) ? AD(1.0) : AD(0.0);
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      AD s = y[i];
      for (std::size_t k = 0; k < i; ++k) {
        s -= l[i][k] * y[k];
      }
      y[i] = s / l[i][i];
    }
    for (int i = static_cast<int>(kSunD) - 1; i >= 0; --i) {
      AD s = y[static_cast<std::size_t>(i)];
      for (std::size_t k = static_cast<std::size_t>(i) + 1; k < kSunD; ++k) {
        s -= l[k][static_cast<std::size_t>(i)] * y[k];
      }
      y[static_cast<std::size_t>(i)] = s / l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      omega_inv[i][col] = y[i];
    }
  }

  ChAD d_oinv = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        s += dp.Delta[k][i] * omega_inv[k][j];
      }
      d_oinv[i][j] = s;
    }
  }

  ChAD lambda_ad = dp.Gamma;
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        s += d_oinv[i][k] * dp.Delta[k][j];
      }
      lambda_ad[i][j] -= s;
    }
  }

  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      lambda[i][j] = Value(lambda_ad[i][j]);
    }
  }
  Mat3 omega{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      omega[i][j] = Value(dp.Omega[i][j]);
    }
  }
  return omega;
}

AD log_mvn3(const std::array<AD, kSunD>& x, const SunDpAD& dp) {
  ChAD l{};
  chol_lower_ad3(dp.Omega, l);
  AD log_det = AD(0.0);
  for (std::size_t i = 0; i < kSunD; ++i) {
    log_det += CppAD::log(l[i][i]);
  }
  log_det *= AD(2.0);

  std::array<AD, kSunD> y{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    y[i] = x[i] - dp.xi[i];
  }
  for (std::size_t i = 0; i < kSunD; ++i) {
    AD s = y[i];
    for (std::size_t k = 0; k < i; ++k) {
      s -= l[i][k] * y[k];
    }
    y[i] = s / l[i][i];
  }
  AD quad = AD(0.0);
  for (std::size_t i = 0; i < kSunD; ++i) {
    quad += y[i] * y[i];
  }
  return AD(-1.5) * AD(kLog2Pi) - AD(0.5) * log_det - AD(0.5) * quad;
}

ChAD lambda_from_dp(const SunDpAD& dp) {
  ChAD omega_bar = cov2cor3(dp.Omega);
  ChAD omega_inv = zeros_ad3();
  ChAD l{};
  chol_lower_ad3(omega_bar, l);
  for (std::size_t col = 0; col < kSunD; ++col) {
    std::array<AD, kSunD> y{};
    for (std::size_t i = 0; i < kSunD; ++i) {
      y[i] = (i == col) ? AD(1.0) : AD(0.0);
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      AD s = y[i];
      for (std::size_t k = 0; k < i; ++k) {
        s -= l[i][k] * y[k];
      }
      y[i] = s / l[i][i];
    }
    for (int i = static_cast<int>(kSunD) - 1; i >= 0; --i) {
      AD s = y[static_cast<std::size_t>(i)];
      for (std::size_t k = static_cast<std::size_t>(i) + 1; k < kSunD; ++k) {
        s -= l[k][static_cast<std::size_t>(i)] * y[k];
      }
      y[static_cast<std::size_t>(i)] = s / l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      omega_inv[i][col] = y[i];
    }
  }

  ChAD d_oinv = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        s += dp.Delta[k][i] * omega_inv[k][j];
      }
      d_oinv[i][j] = s;
    }
  }

  ChAD lambda = dp.Gamma;
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        s += d_oinv[i][k] * dp.Delta[k][j];
      }
      lambda[i][j] -= s;
    }
  }
  return lambda;
}

std::array<AD, kSunD> alpha_from_x(const std::array<AD, kSunD>& x, const SunDpAD& dp) {
  ChAD omega_bar = cov2cor3(dp.Omega);
  ChAD omega_inv = zeros_ad3();
  ChAD l{};
  chol_lower_ad3(omega_bar, l);
  for (std::size_t col = 0; col < kSunD; ++col) {
    std::array<AD, kSunD> y{};
    for (std::size_t i = 0; i < kSunD; ++i) {
      y[i] = (i == col) ? AD(1.0) : AD(0.0);
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      AD s = y[i];
      for (std::size_t k = 0; k < i; ++k) {
        s -= l[i][k] * y[k];
      }
      y[i] = s / l[i][i];
    }
    for (int i = static_cast<int>(kSunD) - 1; i >= 0; --i) {
      AD s = y[static_cast<std::size_t>(i)];
      for (std::size_t k = static_cast<std::size_t>(i) + 1; k < kSunD; ++k) {
        s -= l[k][static_cast<std::size_t>(i)] * y[k];
      }
      y[static_cast<std::size_t>(i)] = s / l[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)];
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      omega_inv[i][col] = y[i];
    }
  }

  std::array<AD, kSunD> z{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    const AD omega_i = CppAD::sqrt(dp.Omega[i][i]);
    z[i] = (x[i] - dp.xi[i]) / omega_i;
  }

  ChAD d_oinv = zeros_ad3();
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      AD s = AD(0.0);
      for (std::size_t k = 0; k < kSunD; ++k) {
        s += dp.Delta[k][i] * omega_inv[k][j];
      }
      d_oinv[i][j] = s;
    }
  }

  std::array<AD, kSunD> alpha{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    AD s = AD(0.0);
    for (std::size_t j = 0; j < kSunD; ++j) {
      s += d_oinv[i][j] * z[j];
    }
    alpha[i] = s;
  }
  return alpha;
}

void pack_genz_ad(
  const ChAD& sigma,
  std::array<AD, kSunD>& scale,
  std::vector<AD>& vech_ch) {

  ChAD corr = cov2cor3(sigma);
  ChAD ch{};
  chol_lower_ad3(corr, ch);
  for (std::size_t i = 0; i < kSunD; ++i) {
    const AD s2 = CppAD::CondExpGt(sigma[i][i], AD(1e-12), sigma[i][i], AD(1e-12));
    scale[i] = CppAD::sqrt(s2);
  }
  vech_ch.assign(vech_size(kSunD), AD(0.0));
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      vech_ch[vech_index(i, j)] = ch[i][j];
    }
  }
}

AD sun_obs_contrib_ad(
  const std::array<AD, kSunD>& x,
  const SunDpAD& dp,
  const ChAD& lambda,
  const std::array<AD, kSunD>& scale_l,
  const std::vector<AD>& vech_l,
  MvnTape& p1_tape) {

  const AD log_pdf = log_mvn3(x, dp);
  const std::array<AD, kSunD> alpha = alpha_from_x(x, dp);
  std::vector<AD> alpha_v(alpha.begin(), alpha.end());
  std::vector<AD> mean_v(kSunD, AD(0.0));
  std::vector<AD> scale_l_v(scale_l.begin(), scale_l.end());
  const AD p1 = detail::pmvn_ad(0, p1_tape, alpha_v, mean_v, scale_l_v, vech_l);
  const AD eps = AD(1e-300);
  return log_pdf + CppAD::log(p1 + eps);
}

AD sun_obs_loglik_ad(
  const std::array<AD, kSunD>& x,
  const std::vector<AD>& par,
  MvnTape& p1_tape) {

  SunDpAD dp = make_sun_dp_ad(par);
  const ChAD lambda = lambda_from_dp(dp);
  std::array<AD, kSunD> scale_l{};
  std::vector<AD> vech_l;
  pack_genz_ad(lambda, scale_l, vech_l);
  return sun_obs_contrib_ad(x, dp, lambda, scale_l, vech_l, p1_tape);
}

AD sun_log_p2_ad(const std::vector<AD>& par, MvnTape& p2_tape) {
  SunDpAD dp = make_sun_dp_ad(par);
  std::array<AD, kSunD> scale_g{};
  std::vector<AD> vech_g;
  pack_genz_ad(dp.Gamma, scale_g, vech_g);

  std::vector<AD> mean_v(kSunD, AD(0.0));
  std::vector<AD> upper_zero(kSunD, AD(0.0));
  std::vector<AD> scale_g_v(scale_g.begin(), scale_g.end());
  const AD p2 = detail::pmvn_ad(1, p2_tape, upper_zero, mean_v, scale_g_v, vech_g);
  const AD eps = AD(1e-300);
  return CppAD::log(p2 + eps);
}

std::vector<std::vector<double>> mat_from_mat3(const Mat3& m) {
  std::vector<std::vector<double>> out(kSunD, std::vector<double>(kSunD));
  for (std::size_t i = 0; i < kSunD; ++i) {
    for (std::size_t j = 0; j < kSunD; ++j) {
      out[i][j] = m[i][j];
    }
  }
  return out;
}

void setup_adfun_sparsity(SunAdfunPack& pack, const std::vector<double>& par_seed) {
  pack.x_seed = detail::to_cppad_vector(par_seed);
  std::vector<std::size_t> par_subset(kSunNPar);
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

void build_obs_shard_fun(SunObsShard& shard) {
  detail::init_qnorm_atomic();
  detail::init_pmvn_atomic();
  detail::set_pmvn_atomic_tapes(&shard.p1_tape, nullptr);

  std::array<AD, kSunD> x_ad{};
  for (std::size_t i = 0; i < kSunD; ++i) {
    x_ad[i] = shard.x[i];
  }

  std::vector<AD> par_ad(kSunNPar);
  for (std::size_t i = 0; i < kSunNPar; ++i) {
    par_ad[i] = AD(0.0);
  }
  CppAD::Independent(par_ad);

  AD y = sun_obs_loglik_ad(x_ad, par_ad, shard.p1_tape);
  std::vector<AD> out(1);
  out[0] = y;
  shard.adfun.fun = CppAD::ADFun<double>(par_ad, out);
}

void build_p2_shard_fun(SunP2Shard& shard) {
  detail::init_qnorm_atomic();
  detail::init_pmvn_atomic();
  detail::set_pmvn_atomic_tapes(nullptr, &shard.p2_tape);

  std::vector<AD> par_ad(kSunNPar);
  for (std::size_t i = 0; i < kSunNPar; ++i) {
    par_ad[i] = AD(0.0);
  }
  CppAD::Independent(par_ad);

  AD y = sun_log_p2_ad(par_ad, shard.p2_tape);
  std::vector<AD> out(1);
  out[0] = y;
  shard.adfun.fun = CppAD::ADFun<double>(par_ad, out);
}

SunResult eval_adfun_log(
  SunAdfunPack& pack,
  MvnTape* p1_tape,
  MvnTape* p2_tape,
  const std::vector<double>& par,
  int deriv) {

  SunResult out;
  out.gradient.assign(kSunNPar, 0.0);
  out.hessian.assign(kSunNPar, std::vector<double>(kSunNPar, 0.0));

  detail::init_qnorm_atomic();
  detail::init_pmvn_atomic();
  detail::set_pmvn_atomic_tapes(p1_tape, p2_tape);

  const CppAD::vector<double> u = detail::to_cppad_vector(par);
  const double raw = pack.fun.Forward(0, u)[0];
  out.value = raw;
  out.error = 0.0;

  if (deriv < 1) {
    return out;
  }

  // Dense reverse: atomic_pmvn does not implement jac_sparsity, so
  // for_jac_sparsity / sparse_jac_rev can drop dependence on parameters that
  // only enter through the orthant CDF (L*, a, b, c), yielding exact zeros.
  CppAD::vector<double> w_rev(1);
  w_rev[0] = 1.0;
  const CppAD::vector<double> dw = pack.fun.Reverse(1, w_rev);
  for (std::size_t j = 0; j < kSunNPar; ++j) {
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

void apply_log_scale(SunResult& out, bool log_scale, int deriv) {
  if (log_scale) {
    return;
  }
  const std::vector<double> grad_log = out.gradient;
  const std::vector<std::vector<double>> h_log = out.hessian;
  const double f = std::exp(out.value);
  out.value = f;
  if (deriv >= 1) {
    for (std::size_t i = 0; i < kSunNPar; ++i) {
      out.gradient[i] = f * grad_log[i];
    }
  }
  if (deriv >= 2) {
    for (std::size_t i = 0; i < kSunNPar; ++i) {
      for (std::size_t j = 0; j < kSunNPar; ++j) {
        out.hessian[i][j] = f * (h_log[i][j] + grad_log[i] * grad_log[j]);
      }
    }
  }
}

void accumulate_result(
  SunResult& total,
  const SunResult& part,
  double scale,
  int deriv) {

  total.value += scale * part.value;
  if (deriv >= 1) {
    for (std::size_t i = 0; i < kSunNPar; ++i) {
      total.gradient[i] += scale * part.gradient[i];
    }
  }
  if (deriv >= 2) {
    for (std::size_t i = 0; i < kSunNPar; ++i) {
      for (std::size_t j = 0; j < kSunNPar; ++j) {
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
  // Cap by hardware concurrency, not omp_get_max_threads(). The latter
  // reflects the last omp_set_num_threads() call; CppadParallelScope teardown
  // sets that to 1, which would permanently serialize later requests.
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

SunTapeBundle create_sun_bundle(
  const std::vector<std::vector<double>>& x_rows,
  const std::vector<double>& par_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  int n_threads,
  const std::vector<double>& weights) {

  if (x_rows.empty() || par_seed.size() != kSunNPar) {
    throw std::runtime_error("invalid seed dimensions for SUN(3,3) tape");
  }
  for (const auto& row : x_rows) {
    if (row.size() != kSunD) {
      throw std::runtime_error("each observation must have length 3");
    }
  }
  if (!weights.empty() && weights.size() != x_rows.size()) {
    throw std::runtime_error("weights must have length equal to number of rows in x");
  }

  SunTapeBundle bundle;
  bundle.n_obs = x_rows.size();
  bundle.n_threads = resolve_n_threads(n_threads, 1);
  if (weights.empty()) {
    bundle.weights.assign(x_rows.size(), 1.0);
  } else {
    bundle.weights = weights;
  }
  bundle.weight_sum = 0.0;
  for (double w : bundle.weights) {
    bundle.weight_sum += w;
  }

  Mat3 lambda{};
  (void)sun_sigma_double(par_seed, lambda);
  Mat3 gamma{};
  {
    std::vector<AD> par_ad(kSunNPar);
    for (std::size_t i = 0; i < kSunNPar; ++i) {
      par_ad[i] = AD(par_seed[i]);
    }
    SunDpAD dp = make_sun_dp_ad(par_ad);
    for (std::size_t i = 0; i < kSunD; ++i) {
      for (std::size_t j = 0; j < kSunD; ++j) {
        gamma[i][j] = Value(dp.Gamma[i][j]);
      }
    }
  }

  std::vector<double> lower(kSunD, kNegInf);
  std::vector<double> mean_zero(kSunD, 0.0);

  bundle.p2.p2_tape = create_mvn_tape(
    lower, mean_zero, mean_zero, mat_from_mat3(gamma),
    n_points, n_shifts, seed + 1U, /*value_only=*/true);
  build_p2_shard_fun(bundle.p2);
  setup_adfun_sparsity(bundle.p2.adfun, par_seed);

  bundle.shards.reserve(x_rows.size());
  for (std::size_t obs = 0; obs < x_rows.size(); ++obs) {
    SunObsShard shard;
    for (std::size_t i = 0; i < kSunD; ++i) {
      shard.x[i] = x_rows[obs][i];
    }

    std::vector<AD> par_ad(kSunNPar);
    std::array<AD, kSunD> x_ad{};
    for (std::size_t i = 0; i < kSunNPar; ++i) {
      par_ad[i] = AD(par_seed[i]);
    }
    for (std::size_t i = 0; i < kSunD; ++i) {
      x_ad[i] = AD(x_rows[obs][i]);
    }
    SunDpAD dp = make_sun_dp_ad(par_ad);
    std::array<AD, kSunD> alpha_ad = alpha_from_x(x_ad, dp);
    std::vector<double> upper_alpha(kSunD);
    for (std::size_t i = 0; i < kSunD; ++i) {
      upper_alpha[i] = Value(alpha_ad[i]);
    }

    shard.p1_tape = create_mvn_tape(
      lower, upper_alpha, mean_zero, mat_from_mat3(lambda),
      n_points, n_shifts, seed + 2U + static_cast<unsigned int>(obs),
      /*value_only=*/true);

    build_obs_shard_fun(shard);
    setup_adfun_sparsity(shard.adfun, par_seed);
    bundle.shards.push_back(std::move(shard));
  }

  return bundle;
}

SunResult eval_sun_bundle(
  SunTapeBundle& bundle,
  const std::vector<double>& par,
  bool log_scale,
  int deriv,
  int n_threads,
  bool manage_parallel_scope) {

  if (par.size() != kSunNPar) {
    throw std::runtime_error("par must have length 21");
  }
  if (deriv < 0 || deriv > 2) {
    throw std::runtime_error("deriv must be 0, 1, or 2");
  }

  n_threads = resolve_n_threads(n_threads, bundle.n_threads);

  SunResult total;
  total.gradient.assign(kSunNPar, 0.0);
  total.hessian.assign(kSunNPar, std::vector<double>(kSunNPar, 0.0));

  // Optional RAII scope; trustOptim driver owns one scope for the whole loop.
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
      shard_grads.assign(n_shards, std::vector<double>(kSunNPar, 0.0));
    }
    if (deriv >= 2) {
      shard_hess.assign(
        n_shards,
        std::vector<std::vector<double>>(kSunNPar, std::vector<double>(kSunNPar, 0.0)));
    }

#pragma omp parallel num_threads(n_threads)
    {
#pragma omp for schedule(static)
      for (std::size_t i = 0; i < n_shards; ++i) {
        SunResult shard_res = eval_adfun_log(
          bundle.shards[i].adfun,
          &bundle.shards[i].p1_tape,
          nullptr,
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
        for (std::size_t j = 0; j < kSunNPar; ++j) {
          total.gradient[j] += wi * shard_grads[i][j];
        }
      }
      if (deriv >= 2) {
        for (std::size_t r = 0; r < kSunNPar; ++r) {
          for (std::size_t c = 0; c < kSunNPar; ++c) {
            total.hessian[r][c] += wi * shard_hess[i][r][c];
          }
        }
      }
    }
  } else
#endif
  {
    for (std::size_t i = 0; i < n_shards; ++i) {
      SunResult shard_res = eval_adfun_log(
        bundle.shards[i].adfun,
        &bundle.shards[i].p1_tape,
        nullptr,
        par,
        deriv);
      accumulate_result(total, shard_res, bundle.weights[i], deriv);
    }
  }

  SunResult p2_res = eval_adfun_log(
    bundle.p2.adfun,
    nullptr,
    &bundle.p2.p2_tape,
    par,
    deriv);
  accumulate_result(total, p2_res, -bundle.weight_sum, deriv);

  apply_log_scale(total, log_scale, deriv);
  return total;
}

}  // namespace admvn
