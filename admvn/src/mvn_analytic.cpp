#include "mvn_analytic.hpp"

#include <Rcpp.h>
#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace admvn {
namespace {

constexpr double kTwoPi = 2.0 * 3.14159265358979323846;

inline bool is_neg_inf(double x) {
  // Accept true -Inf and the historical Genz sentinel (~-1e50).
  return (!std::isfinite(x) && x < 0.0) || x <= -1e20;
}

inline double clamp_cor(double r) {
  if (r > 0.999999) {
    return 0.999999;
  }
  if (r < -0.999999) {
    return -0.999999;
  }
  return r;
}

inline double dnorm0(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

inline double pnorm0(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
}

// 16-point Gauss-Legendre nodes/weights on (-1, 1).
const double kGl16X[8] = {
  0.095012509837637440185, 0.281603550779258913230, 0.458016777657227386342,
  0.617876244402643748447, 0.755404408355003033895, 0.865631202387831743880,
  0.944575023073232576078, 0.989400934991649932596
};
const double kGl16W[8] = {
  0.189450610455068496285, 0.182603415044923588867, 0.169156519395002538189,
  0.149595988816576732081, 0.124628971255533872052, 0.095158511682492784810,
  0.062253523938647892863, 0.027152459411754094852
};

inline double clamp01(double p) {
  if (p < 0.0) {
    return 0.0;
  }
  if (p > 1.0) {
    return 1.0;
  }
  return p;
}

inline double finite_or_clip(double x) {
  if (std::isfinite(x)) {
    return x;
  }
  return x > 0.0 ? 8.0 : -8.0;
}

// Phi_2(h,k;r) = int_{-inf}^h phi(x) Phi((k-r x)/sqrt(1-r^2)) dx
// (Drezner/Genz-style Gauss-Legendre; sn/mnormt bivariate path).
double pbvnorm(double h, double k, double r) {
  r = clamp_cor(r);
  if ((!std::isfinite(h) && h < 0.0) || (!std::isfinite(k) && k < 0.0)) {
    return 0.0;
  }
  if (!std::isfinite(h) && h > 0.0 && !std::isfinite(k) && k > 0.0) {
    return 1.0;
  }
  if (std::abs(r) < 1e-15) {
    return pnorm0(h) * pnorm0(k);
  }
  if (r > 0.99999) {
    return pnorm0(std::min(h, k));
  }
  if (r < -0.99999) {
    return std::max(0.0, pnorm0(h) + pnorm0(k) - 1.0);
  }

  const double hh = finite_or_clip(h);
  const double kk = finite_or_clip(k);
  const double s = std::sqrt(std::max(1.0 - r * r, 1e-16));

  const double lo = std::min(-8.0, hh - 8.0);
  const double mid = 0.5 * (hh + lo);
  const double half = 0.5 * (hh - lo);
  double acc = 0.0;
  for (int i = 0; i < 8; ++i) {
    const double x1 = mid + half * kGl16X[i];
    const double x2 = mid - half * kGl16X[i];
    acc += kGl16W[i] * (
      dnorm0(x1) * pnorm0((kk - r * x1) / s) +
      dnorm0(x2) * pnorm0((kk - r * x2) / s));
  }
  return clamp01(acc * half);
}

// Phi_3(h; R) via int_{-inf}^{h1} phi(x) Phi_2(h2|x, h3|x; r23|1) dx
// (same reduction used by specialized trivariate normal CDF routines).
double ptvnorm(
  double h1, double h2, double h3,
  double r12, double r13, double r23) {

  if ((!std::isfinite(h1) && h1 < 0.0) ||
      (!std::isfinite(h2) && h2 < 0.0) ||
      (!std::isfinite(h3) && h3 < 0.0)) {
    return 0.0;
  }
  if (!std::isfinite(h1) && h1 > 0.0 &&
      !std::isfinite(h2) && h2 > 0.0 &&
      !std::isfinite(h3) && h3 > 0.0) {
    return 1.0;
  }

  r12 = clamp_cor(r12);
  r13 = clamp_cor(r13);
  r23 = clamp_cor(r23);

  // Degenerate / nearly singular pairwise correlations.
  if (std::abs(r12) > 0.99999 && std::abs(r13) > 0.99999 && std::abs(r23) > 0.99999) {
    if (r12 > 0.0 && r13 > 0.0 && r23 > 0.0) {
      return pnorm0(std::min(h1, std::min(h2, h3)));
    }
  }

  const double hh = finite_or_clip(h1);
  const double lo = std::min(-8.0, hh - 8.0);
  const double mid = 0.5 * (hh + lo);
  const double half = 0.5 * (hh - lo);

  auto integrand = [&](double x) {
    const double v2 = std::max(1.0 - r12 * r12, 1e-16);
    const double v3 = std::max(1.0 - r13 * r13, 1e-16);
    const double s2 = std::sqrt(v2);
    const double s3 = std::sqrt(v3);
    const double a2 = (h2 - r12 * x) / s2;
    const double a3 = (h3 - r13 * x) / s3;
    const double cov = r23 - r12 * r13;
    const double rho = clamp_cor(cov / (s2 * s3));
    return dnorm0(x) * pbvnorm(a2, a3, rho);
  };

  double acc = 0.0;
  for (int i = 0; i < 8; ++i) {
    acc += kGl16W[i] * (integrand(mid + half * kGl16X[i]) +
                        integrand(mid - half * kGl16X[i]));
  }
  return clamp01(acc * half);
}

double dbvnorm(double x, double y, double r) {
  r = clamp_cor(r);
  const double s = 1.0 - r * r;
  const double q = (x * x - 2.0 * r * x * y + y * y) / s;
  return std::exp(-0.5 * q) / (kTwoPi * std::sqrt(s));
}

bool all_lower_neg_inf(const std::vector<double>& lower, std::size_t n) {
  if (lower.size() != n) {
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!is_neg_inf(lower[i])) {
      return false;
    }
  }
  return true;
}

// Conditionally evaluate P(X_{-k} <= upper_{-k} | X_k = upper_k)
// for X ~ N(mean, Sigma).
double cond_cdf_drop_one_cov(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  std::size_t k) {

  const std::size_t n = upper.size();
  if (n == 1) {
    return 1.0;
  }

  std::vector<std::size_t> idx;
  idx.reserve(n - 1);
  for (std::size_t i = 0; i < n; ++i) {
    if (i != k) {
      idx.push_back(i);
    }
  }
  const std::size_t m = idx.size();
  const double skk = sigma[k][k];
  if (!(skk > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  std::vector<double> cmean(m);
  std::vector<std::vector<double>> ccov(m, std::vector<double>(m));
  for (std::size_t ii = 0; ii < m; ++ii) {
    cmean[ii] = mean[idx[ii]] + sigma[idx[ii]][k] / skk * (upper[k] - mean[k]);
    for (std::size_t jj = 0; jj < m; ++jj) {
      ccov[ii][jj] =
        sigma[idx[ii]][idx[jj]] - sigma[idx[ii]][k] * sigma[idx[jj]][k] / skk;
    }
  }

  if (m == 1) {
    const double sd = std::sqrt(std::max(ccov[0][0], 1e-16));
    return pnorm0((upper[idx[0]] - cmean[0]) / sd);
  }

  if (m == 2) {
    const double s0 = std::sqrt(std::max(ccov[0][0], 1e-16));
    const double s1 = std::sqrt(std::max(ccov[1][1], 1e-16));
    const double rho = clamp_cor(ccov[0][1] / (s0 * s1));
    return pbvnorm(
      (upper[idx[0]] - cmean[0]) / s0,
      (upper[idx[1]] - cmean[1]) / s1,
      rho);
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double cond_cdf_drop_two_std(
  const std::vector<double>& a,
  const std::vector<std::vector<double>>& R,
  std::size_t i,
  std::size_t j) {

  const std::size_t n = a.size();
  if (n == 2) {
    return 1.0;
  }
  std::size_t k = n;
  for (std::size_t t = 0; t < n; ++t) {
    if (t != i && t != j) {
      k = t;
      break;
    }
  }
  if (k >= n) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double rho = clamp_cor(R[i][j]);
  const double det = std::max(1.0 - rho * rho, 1e-16);
  const double pii = 1.0 / det;
  const double pij = -rho / det;
  const double pjj = 1.0 / det;

  const double mu =
    R[k][i] * (pii * a[i] + pij * a[j]) +
    R[k][j] * (pij * a[i] + pjj * a[j]);
  const double var =
    R[k][k] -
    (R[k][i] * (pii * R[k][i] + pij * R[k][j]) +
     R[k][j] * (pij * R[k][i] + pjj * R[k][j]));
  return pnorm0((a[k] - mu) / std::sqrt(std::max(var, 1e-16)));
}

std::vector<std::vector<double>> sigma_from_genz(const GenzPack& genz) {
  const std::size_t n = genz.scale.size();
  std::vector<std::vector<double>> sigma(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double r = 0.0;
      for (std::size_t k = 0; k <= j; ++k) {
        r += genz.ch[i][k] * genz.ch[j][k];
      }
      sigma[i][j] = genz.scale[i] * r * genz.scale[j];
      sigma[j][i] = sigma[i][j];
    }
  }
  return sigma;
}

std::vector<std::vector<double>> corr_from_sigma(
  const std::vector<std::vector<double>>& sigma,
  std::vector<double>& scale) {

  const std::size_t n = sigma.size();
  scale.resize(n);
  std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    scale[i] = std::sqrt(std::max(sigma[i][i], 1e-16));
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      R[i][j] = sigma[i][j] / (scale[i] * scale[j]);
    }
  }
  return R;
}

bool fill_upper_mean_grad(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  std::vector<double>& grad_upper,
  std::vector<double>& grad_mean) {

  const std::size_t n = upper.size();
  grad_upper.assign(n, 0.0);
  grad_mean.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double sd = std::sqrt(std::max(sigma[i][i], 1e-16));
    const double cond = cond_cdf_drop_one_cov(upper, mean, sigma, i);
    if (!std::isfinite(cond)) {
      return false;
    }
    const double du = R::dnorm(upper[i], mean[i], sd, 0) * cond;
    grad_upper[i] = du;
    grad_mean[i] = -du;
  }
  return true;
}

}  // namespace

double pmvn_cdf_std(
  const std::vector<double>& h,
  const std::vector<std::vector<double>>& R) {

  const std::size_t n = h.size();
  if (n == 0 || n > 3 || R.size() != n) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (R[i].size() != n) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  if (n == 1) {
    return pnorm0(h[0]);
  }
  if (n == 2) {
    return pbvnorm(h[0], h[1], R[0][1]);
  }
  // n == 3
  return ptvnorm(h[0], h[1], h[2], R[0][1], R[0][2], R[1][2]);
}

double pmvn_cdf_special(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  const std::vector<double>& lower) {

  const std::size_t n = upper.size();
  if (n == 0 || n > 3 || mean.size() != n || sigma.size() != n ||
      !all_lower_neg_inf(lower, n)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (sigma[i].size() != n || !(sigma[i][i] > 0.0)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  std::vector<double> scale;
  const std::vector<std::vector<double>> R = corr_from_sigma(sigma, scale);
  std::vector<double> h(n);
  for (std::size_t i = 0; i < n; ++i) {
    h[i] = (upper[i] - mean[i]) / scale[i];
  }
  return pmvn_cdf_std(h, R);
}

bool analytic_mvn_upper_mean_grad(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  const std::vector<double>& lower,
  std::vector<double>& grad_upper,
  std::vector<double>& grad_mean) {

  const std::size_t n = upper.size();
  if (n == 0 || n > 3 || mean.size() != n || sigma.size() != n ||
      !all_lower_neg_inf(lower, n)) {
    return false;
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (sigma[i].size() != n) {
      return false;
    }
  }
  return fill_upper_mean_grad(upper, mean, sigma, grad_upper, grad_mean);
}

std::vector<double> analytic_mvn_domain_grad(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz,
  const std::vector<double>& lower) {

  const std::size_t n = upper.size();
  if (n == 0 || n > 3 || mean.size() != n || genz.scale.size() != n ||
      genz.ch.size() != n || !all_lower_neg_inf(lower, n)) {
    return {};
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!(genz.scale[i] > 0.0) || genz.ch[i].size() != n) {
      return {};
    }
  }

  const std::vector<std::vector<double>> sigma = sigma_from_genz(genz);

  std::vector<double> grad_upper;
  std::vector<double> grad_mean;
  if (!fill_upper_mean_grad(upper, mean, sigma, grad_upper, grad_mean)) {
    return {};
  }

  // Standardized Plackett derivatives for dF/d rho
  std::vector<double> scale;
  const std::vector<std::vector<double>> R = corr_from_sigma(sigma, scale);
  std::vector<double> a(n);
  for (std::size_t i = 0; i < n; ++i) {
    a[i] = (upper[i] - mean[i]) / scale[i];
  }

  std::vector<std::vector<double>> dphi_dr(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const double rho = clamp_cor(R[i][j]);
      const double cond = cond_cdf_drop_two_std(a, R, i, j);
      if (!std::isfinite(cond)) {
        return {};
      }
      dphi_dr[i][j] = dbvnorm(a[i], a[j], rho) * cond;
      dphi_dr[j][i] = dphi_dr[i][j];
    }
  }

  // dPhi/da from upper grads: dF/du_i = dPhi/da_i / scale_i
  std::vector<double> dphi_da(n);
  for (std::size_t i = 0; i < n; ++i) {
    dphi_da[i] = grad_upper[i] * scale[i];
  }

  const std::size_t n_dom = domain_size(n);
  std::vector<double> grad(n_dom, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    grad[i] = grad_upper[i];
    grad[n + i] = grad_mean[i];
    // scale in domain is genz.scale; a_i = (u-mu)/genz.scale
    grad[2 * n + i] = dphi_da[i] * (-a[i] / genz.scale[i]);
  }

  const std::size_t off_ch = 3 * n;
  for (std::size_t a_idx = 0; a_idx < n; ++a_idx) {
    for (std::size_t b = 0; b <= a_idx; ++b) {
      double g = 0.0;
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
          double dr = 0.0;
          if (a_idx == i && b <= j) {
            dr += genz.ch[j][b];
          }
          if (a_idx == j && b <= i) {
            dr += genz.ch[i][b];
          }
          g += dphi_dr[i][j] * dr;
        }
      }
      grad[off_ch + vech_index(a_idx, b)] = g;
    }
  }

  return grad;
}

}  // namespace admvn
