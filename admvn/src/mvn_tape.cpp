#include "mvn_tape.hpp"

#include "mvn_analytic.hpp"
#include "mvn_sparsity.hpp"
#include "qnorm_atomic.hpp"

#include <Rcpp.h>
#include <Rmath.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

namespace admvn {
namespace {

constexpr double kNine = 9.0;
constexpr double kEp = 1e-10;

inline bool is_neg_inf(double x) {
  return !std::isfinite(x) && x < 0.0;
}

inline bool is_pos_inf(double x) {
  return !std::isfinite(x) && x > 0.0;
}

inline double cdf_std(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
}

inline double cdf_scaled(double x, double ct) {
  if (x < -kNine * ct) {
    return 0.0;
  }
  if (x > kNine * ct) {
    return 1.0;
  }
  return cdf_std(x / ct);
}

std::vector<double> primes_sqrt(std::size_t n) {
  std::vector<double> out;
  if (n <= 1) {
    return out;
  }
  out.reserve(n - 1);
  std::size_t count = 0;
  int candidate = 2;
  while (count < n - 1) {
    bool prime = true;
    for (int p = 2; p * p <= candidate; ++p) {
      if (candidate % p == 0) {
        prime = false;
        break;
      }
    }
    if (prime) {
      if (count > 0) {
        out.push_back(std::sqrt(static_cast<double>(candidate)));
      }
      ++count;
    }
    ++candidate;
  }
  return out;
}

void chlrdr(
  std::vector<std::vector<double>>& c,
  std::vector<double>& ap,
  std::vector<double>& bp,
  std::vector<int>& perm) {

  const std::size_t n = c.size();
  const double eps = std::numeric_limits<double>::epsilon();
  std::vector<double> d(n);

  for (std::size_t i = 0; i < n; ++i) {
    d[i] = std::sqrt(c[i][i]);
  }

  for (std::size_t i = 0; i < n; ++i) {
    if (d[i] > 0.0) {
      const double di = d[i];
      for (std::size_t j = 0; j < n; ++j) {
        c[i][j] /= di;
        c[j][i] /= di;
      }
      ap[i] /= di;
      bp[i] /= di;
    }
  }

  std::vector<double> y(n, 0.0);

  for (std::size_t k = 0; k < n; ++k) {
    std::size_t im = k;
    double ckk = 0.0;
    double dem = 1.0;
    double am = 0.0;
    double bm = 0.0;

    for (std::size_t i = k; i < n; ++i) {
      if (c[i][i] > eps) {
        const double cii = std::sqrt(std::max(c[i][i], 0.0));
        double s = 0.0;
        if (i > 0 && k > 0) {
          for (std::size_t j = 0; j < k; ++j) {
            s += c[i][j] * y[j];
          }
        }
        const double ai = (ap[i] - s) / cii;
        const double bi = (bp[i] - s) / cii;
        const double de = cdf_std(bi) - cdf_std(ai);
        if (de <= dem) {
          ckk = cii;
          dem = de;
          am = ai;
          bm = bi;
          im = i;
        }
      }
    }

    if (im > k) {
      c[im][im] = c[k][k];
      std::swap(ap[im], ap[k]);
      std::swap(bp[im], bp[k]);
      std::swap(perm[im], perm[k]);
      if (k > 0) {
        for (std::size_t j = 0; j < k; ++j) {
          std::swap(c[im][j], c[k][j]);
        }
      }
      if (im < n - 1) {
        for (std::size_t j = im + 1; j < n; ++j) {
          std::swap(c[j][im], c[j][k]);
        }
      }
      if (k <= n - 2 && im <= n - 1) {
        for (std::size_t j = k + 1; j < im; ++j) {
          std::swap(c[j][k], c[im][j]);
        }
      }
    }

    if (k < n - 1) {
      for (std::size_t j = k + 1; j < n; ++j) {
        c[k][j] = 0.0;
      }
    }

    if (ckk > static_cast<double>(k + 1) * kEp) {
      c[k][k] = ckk;
      for (std::size_t i = k + 1; i < n; ++i) {
        c[i][k] /= ckk;
        for (std::size_t j = k + 1; j <= i; ++j) {
          c[i][j] -= c[i][k] * c[j][k];
        }
      }
      if (std::abs(dem) > kEp) {
        y[k] = (std::exp(-0.5 * am * am) - std::exp(-0.5 * bm * bm)) /
          (R::dnorm(0.0, 0.0, 1.0, 0) * 2.0 * dem);
      } else if (am < -10.0) {
        y[k] = bm;
      } else if (bm > 10.0) {
        y[k] = am;
      } else {
        y[k] = 0.5 * (am + bm);
      }
    } else {
      for (std::size_t i = k; i < n; ++i) {
        c[i][k] = 0.0;
      }
      y[k] = 0.0;
    }
  }
}

bool chol_lower(
  const std::vector<std::vector<double>>& a,
  std::vector<std::vector<double>>& l) {

  const std::size_t n = a.size();
  l.assign(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      double s = a[i][j];
      for (std::size_t k = 0; k < j; ++k) {
        s -= l[i][k] * l[j][k];
      }
      if (i == j) {
        if (s <= 0.0) {
          return false;
        }
        l[i][j] = std::sqrt(s);
      } else {
        if (l[j][j] == 0.0) {
          return false;
        }
        l[i][j] = s / l[j][j];
      }
    }
  }
  return true;
}

template <typename ADType>
using ChMatrix = std::vector<std::vector<ADType>>;

template <typename ADType>
ADType qmc_integrand(
  const std::vector<ADType>& bs,
  const std::vector<ADType>& as,
  const ChMatrix<ADType>& ch,
  const std::vector<double>& w) {

  const std::size_t n = bs.size();
  const ADType ct0 = ch[0][0];
  const ADType ai0 = as[0];
  const ADType bi0 = bs[0];

  ADType c1 = ADType(0.0);
  if (ai0 > ADType(-kNine) * ct0) {
    if (ai0 < ADType(kNine) * ct0) {
      c1 = detail::pnorm_ad(ai0 / ct0);
    } else {
      c1 = ADType(1.0);
    }
  }

  ADType d1 = ADType(0.0);
  if (bi0 > ADType(-kNine) * ct0) {
    if (bi0 < ADType(kNine) * ct0) {
      d1 = detail::pnorm_ad(bi0 / ct0);
    } else {
      d1 = ADType(1.0);
    }
  }

  ADType c_val = c1;
  ADType dc_val = d1 - c1;
  ADType pv = dc_val;
  std::vector<ADType> y(n - 1, ADType(0.0));

  for (std::size_t i = 1; i < n; ++i) {
    const ADType x = ADType(w[i - 1]);
    y[i - 1] = detail::qnorm_ad(c_val + x * dc_val);

    ADType s = ADType(0.0);
    for (std::size_t k = 0; k < i; ++k) {
      s += ch[i][k] * y[k];
    }

    const ADType ct = ch[i][i];
    const ADType aicnt = as[i] - s;
    const ADType bicnt = bs[i] - s;

    ADType ci = ADType(1.0);
    ADType di = ADType(1.0);

    if (aicnt < ADType(-kNine) * ct) {
      ci = ADType(0.0);
    } else if (CppAD::abs(aicnt) < ADType(kNine) * ct) {
      ci = detail::pnorm_ad(aicnt / ct);
    }

    if (bicnt < ADType(-kNine) * ct) {
      di = ADType(0.0);
    } else if (CppAD::abs(bicnt) < ADType(kNine) * ct) {
      di = detail::pnorm_ad(bicnt / ct);
    }

    dc_val = di - ci;
    pv *= dc_val;
    c_val = ci;
  }

  return pv;
}

double qmc_integrand_double(
  const std::vector<double>& bs,
  const std::vector<double>& as,
  const std::vector<std::vector<double>>& ch,
  const std::vector<double>& w) {

  const std::size_t n = bs.size();
  const double ct0 = ch[0][0];
  double c1 = cdf_scaled(as[0], ct0);
  double d1 = cdf_scaled(bs[0], ct0);
  double c_val = c1;
  double dc_val = d1 - c1;
  double pv = dc_val;
  std::vector<double> y(n - 1, 0.0);

  for (std::size_t i = 1; i < n; ++i) {
    const double x = w[i - 1];
    y[i - 1] = R::qnorm(c_val + x * dc_val, 0.0, 1.0, 1, 0);
    double s = 0.0;
    for (std::size_t k = 0; k < i; ++k) {
      s += ch[i][k] * y[k];
    }
    const double ct = ch[i][i];
    const double aicnt = as[i] - s;
    const double bicnt = bs[i] - s;
    double ci = 1.0;
    double di = 1.0;
    if (aicnt < -kNine * ct) {
      ci = 0.0;
    } else if (std::abs(aicnt) < kNine * ct) {
      ci = cdf_std(aicnt / ct);
    }
    if (bicnt < -kNine * ct) {
      di = 0.0;
    } else if (std::abs(bicnt) < kNine * ct) {
      di = cdf_std(bicnt / ct);
    }
    dc_val = di - ci;
    pv *= dc_val;
    c_val = ci;
  }
  return pv;
}

std::vector<double> make_bs(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<int>& perm,
  const std::vector<double>& scale) {

  const std::size_t n = perm.size();
  std::vector<double> bs(n);
  for (std::size_t i = 0; i < n; ++i) {
    const int j = perm[i];
    bs[i] = (upper[j] - mean[j]) / scale[j];
  }
  return bs;
}

std::vector<double> make_as(
  const std::vector<double>& lower,
  const std::vector<double>& mean,
  const std::vector<int>& perm,
  const std::vector<double>& scale) {

  const std::size_t n = perm.size();
  std::vector<double> as(n);
  for (std::size_t i = 0; i < n; ++i) {
    const int j = perm[i];
    if (is_neg_inf(lower[j])) {
      as[i] = -std::numeric_limits<double>::infinity();
    } else {
      as[i] = (lower[j] - mean[j]) / scale[j];
    }
  }
  return as;
}

ChMatrix<CppAD::AD<double>> unpack_ch_ad(
  const std::vector<CppAD::AD<double>>& x,
  std::size_t n,
  std::size_t ch_offset) {

  ChMatrix<CppAD::AD<double>> ch(n, std::vector<CppAD::AD<double>>(n, CppAD::AD<double>(0.0)));
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      ch[i][j] = x[ch_offset + vech_index(i, j)];
    }
  }
  return ch;
}

void build_taped_fun(
  MvnTape& tape,
  const std::vector<double>& lower,
  const std::vector<double>& domain_seed) {

  detail::init_qnorm_atomic();

  const std::size_t n = tape.n;
  const std::size_t n_dom = tape.n_domain;
  if (domain_seed.size() != n_dom) {
    throw std::runtime_error("domain_seed has wrong length");
  }

  std::vector<CppAD::AD<double>> x_ad(n_dom);
  for (std::size_t i = 0; i < n_dom; ++i) {
    x_ad[i] = domain_seed[i];
  }
  CppAD::Independent(x_ad);

  const std::size_t off_mean = n;
  const std::size_t off_scale = 2 * n;
  const std::size_t off_ch = 3 * n;

  ChMatrix<CppAD::AD<double>> ch_ad = unpack_ch_ad(x_ad, n, off_ch);

  std::vector<CppAD::AD<double>> bs(n);
  std::vector<CppAD::AD<double>> as(n);
  for (std::size_t i = 0; i < n; ++i) {
    const int j = tape.perm[i];
    const CppAD::AD<double> scale_j = x_ad[off_scale + j];
    bs[i] = (x_ad[j] - x_ad[off_mean + j]) / scale_j;
    if (is_neg_inf(lower[j])) {
      as[i] = CppAD::AD<double>(-1e50);
    } else {
      as[i] = (CppAD::AD<double>(lower[j]) - x_ad[off_mean + j]) / scale_j;
    }
  }

  CppAD::AD<double> total = CppAD::AD<double>(0.0);
  const double norm = static_cast<double>(tape.n_shifts * tape.n_points);

  for (std::size_t shift = 0; shift < tape.n_shifts; ++shift) {
    for (std::size_t pt = 0; pt < tape.n_points; ++pt) {
      total += qmc_integrand(bs, as, ch_ad, tape.qmc_w[shift][pt]);
    }
  }
  total /= norm;

  std::vector<CppAD::AD<double>> y(1);
  y[0] = total;
  tape.fun = CppAD::ADFun<double>(x_ad, y);
}

GenzPack pack_genz_ch_sigma_p(
  const std::vector<std::vector<double>>& sigma,
  const std::vector<int>& perm) {

  const std::size_t n = sigma.size();
  GenzPack out;
  out.scale.resize(n);

  std::vector<std::vector<double>> sigma_p(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    out.scale[i] = std::sqrt(sigma[i][i]);
    for (std::size_t j = 0; j < n; ++j) {
      sigma_p[i][j] = sigma[perm[i]][perm[j]];
    }
  }

  std::vector<double> dp(n);
  std::vector<std::vector<double>> corr_p(n, std::vector<double>(n, 0.0));
  for (std::size_t i = 0; i < n; ++i) {
    dp[i] = std::sqrt(sigma_p[i][i]);
    for (std::size_t j = 0; j < n; ++j) {
      corr_p[i][j] = sigma_p[i][j] / (dp[i] * dp[j]);
    }
  }

  if (!chol_lower(corr_p, out.ch)) {
    throw std::runtime_error("pack_genz_ch: Cholesky failed");
  }
  return out;
}

}  // namespace

double eval_mvn_value_double(
  const MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz,
  double* error_out) {

  const std::vector<double> bs = make_bs(upper, mean, tape.perm, genz.scale);
  const std::vector<double> as = make_as(tape.lower, mean, tape.perm, genz.scale);
  double p = 0.0;
  double e = 0.0;

  for (std::size_t shift = 0; shift < tape.n_shifts; ++shift) {
    double sm = 0.0;
    for (std::size_t pt = 0; pt < tape.n_points; ++pt) {
      sm += qmc_integrand_double(bs, as, genz.ch, tape.qmc_w[shift][pt]);
    }
    sm /= static_cast<double>(tape.n_points);
    const double dm = (sm - p) / static_cast<double>(shift + 1);
    p += dm;
    e = (static_cast<double>(shift) * e / static_cast<double>(shift + 1)) + dm * dm;
  }
  if (error_out != nullptr) {
    *error_out = 3.0 * std::sqrt(e);
  }
  return p;
}

std::size_t vech_size(std::size_t n) {
  return n * (n + 1) / 2;
}

std::size_t vech_index(std::size_t i, std::size_t j) {
  return i * (i + 1) / 2 + j;
}

std::size_t domain_size(std::size_t n) {
  return 3 * n + vech_size(n);
}

GenzPack pack_genz_ch(
  const std::vector<std::vector<double>>& sigma,
  const std::vector<int>& perm) {
  return pack_genz_ch_sigma_p(sigma, perm);
}

std::vector<double> pack_domain(
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz) {

  const std::size_t n = upper.size();
  const std::size_t n_dom = domain_size(n);
  std::vector<double> x(n_dom);
  for (std::size_t i = 0; i < n; ++i) {
    x[i] = upper[i];
    x[n + i] = mean[i];
    x[2 * n + i] = genz.scale[i];
  }
  const std::size_t off_ch = 3 * n;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      x[off_ch + vech_index(i, j)] = genz.ch[i][j];
    }
  }
  return x;
}

MvnTape create_mvn_tape(
  const std::vector<double>& lower,
  const std::vector<double>& upper_seed,
  const std::vector<double>& mean_seed,
  const std::vector<std::vector<double>>& sigma_seed,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed,
  bool value_only) {

  const std::size_t n = sigma_seed.size();
  if (n == 0) {
    throw std::runtime_error("sigma must be at least 1x1");
  }

  MvnTape tape;
  tape.n = n;
  tape.n_vech = vech_size(n);
  tape.n_domain = domain_size(n);
  tape.n_points = n_points;
  tape.n_shifts = n_shifts;
  tape.lower = lower;
  tape.value_only = value_only;
  tape.perm.resize(n);
  std::iota(tape.perm.begin(), tape.perm.end(), 0);

  std::vector<std::vector<double>> c = sigma_seed;
  std::vector<double> ap(n);
  std::vector<double> bp(n);
  std::vector<double> scale(n);
  for (std::size_t i = 0; i < n; ++i) {
    scale[i] = std::sqrt(sigma_seed[i][i]);
    ap[i] = is_neg_inf(lower[i]) ? -std::numeric_limits<double>::infinity()
                                 : (lower[i] - mean_seed[i]) / scale[i];
    bp[i] = is_pos_inf(upper_seed[i]) ? std::numeric_limits<double>::infinity()
                                      : (upper_seed[i] - mean_seed[i]) / scale[i];
  }

  chlrdr(c, ap, bp, tape.perm);

  const std::vector<double> q = primes_sqrt(n);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  tape.qmc_w.assign(
    n_shifts,
    std::vector<std::vector<double>>(
      n_points,
      std::vector<double>(n > 1 ? n - 1 : 0, 0.0)));

  for (std::size_t shift = 0; shift < n_shifts; ++shift) {
    const double xr = unif(rng);
    for (std::size_t pt = 0; pt < n_points; ++pt) {
      const double cnt = static_cast<double>(pt + 1);
      for (std::size_t dim = 0; dim + 1 < n; ++dim) {
        double xv = std::fmod(cnt * q[dim] + xr, 1.0);
        if (xv < 0.0) {
          xv += 1.0;
        }
        tape.qmc_w[shift][pt][dim] = std::abs(2.0 * xv - 1.0);
      }
    }
  }

  GenzPack genz;
  genz.scale.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    genz.scale[i] = std::sqrt(sigma_seed[i][i]);
  }
  genz.ch = c;
  const std::vector<double> domain_seed =
    pack_domain(upper_seed, mean_seed, genz);

  tape.x_seed = detail::to_cppad_vector(domain_seed);

  if (!value_only) {
    build_taped_fun(tape, lower, domain_seed);

    std::vector<std::size_t> inner_subset(n);
    std::iota(inner_subset.begin(), inner_subset.end(), 0);

    detail::setup_mvn_sparsity(
      tape.fun,
      tape.x_seed,
      inner_subset,
      tape.pattern_grad,
      tape.pattern_grad_inner,
      tape.pattern_hessian,
      tape.pattern_hessian_inner,
      tape.work_grad,
      tape.work_inner_grad,
      tape.work_hess,
      tape.work_inner_hess,
      tape.unused_pattern,
      tape.w);
  }

  return tape;
}

MvnResult eval_mvn_tape(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  bool inner) {

  const std::size_t n = tape.n;
  MvnResult out;
  out.gradient.assign(n, 0.0);
  out.hessian.assign(n, std::vector<double>(n, 0.0));
  out.gradient_mean.assign(n, 0.0);

  const GenzPack genz = pack_genz_ch(sigma, tape.perm);
  out.value = eval_mvn_value_double(tape, upper, mean, genz, &out.error);

  // Prefer analytic ∂F/∂upper (and mean) from Sigma in original order.
  std::vector<double> g_upper;
  std::vector<double> g_mean;
  if (analytic_mvn_upper_mean_grad(upper, mean, sigma, tape.lower, g_upper, g_mean)) {
    out.gradient = std::move(g_upper);
    out.gradient_mean = std::move(g_mean);
  } else if (tape.value_only) {
    // No AD reverse tape; signal a bad domain (e.g. non-PD Sigma) with NaNs
    // instead of throwing (unsafe under OpenMP).
    const double nan = std::numeric_limits<double>::quiet_NaN();
    out.gradient.assign(n, nan);
    out.gradient_mean.assign(n, nan);
  } else {
    detail::init_qnorm_atomic();
    const CppAD::vector<double> x = detail::to_cppad_vector(pack_domain(upper, mean, genz));
    auto& pattern_grad = inner ? tape.pattern_grad_inner : tape.pattern_grad;
    auto& work_grad = inner ? tape.work_inner_grad : tape.work_grad;
    CppAD::sparse_rc<CppAD::vector<size_t>> jac_pattern;
    tape.fun.sparse_jac_rev(x, pattern_grad, jac_pattern, detail::kJacColor, work_grad);

    if (inner) {
      const auto& cols = pattern_grad.col();
      const auto& vals = pattern_grad.val();
      for (size_t k = 0; k < pattern_grad.nnz(); ++k) {
        out.gradient[cols[k]] += vals[k];
      }
    } else {
      std::vector<double> grad_full(tape.n_domain, 0.0);
      const auto& cols = pattern_grad.col();
      const auto& vals = pattern_grad.val();
      for (size_t k = 0; k < pattern_grad.nnz(); ++k) {
        grad_full[cols[k]] += vals[k];
      }
      for (std::size_t i = 0; i < n; ++i) {
        out.gradient[i] = grad_full[i];
        out.gradient_mean[i] = grad_full[n + i];
      }
    }
  }

  if (!tape.value_only) {
    detail::init_qnorm_atomic();
    const CppAD::vector<double> x = detail::to_cppad_vector(pack_domain(upper, mean, genz));
    auto& pattern_hes = inner ? tape.pattern_hessian_inner : tape.pattern_hessian;
    auto& work_hes = inner ? tape.work_inner_hess : tape.work_hess;
    CppAD::sparse_rc<CppAD::vector<size_t>> hes_pattern;
    tape.fun.sparse_hes(x, tape.w, pattern_hes, hes_pattern, detail::kHessColor, work_hes);

    const auto& rows = pattern_hes.row();
    const auto& cols = pattern_hes.col();
    const auto& vals = pattern_hes.val();
    for (size_t k = 0; k < pattern_hes.nnz(); ++k) {
      const std::size_t i = rows[k];
      const std::size_t j = cols[k];
      if (i < n && j < n) {
        out.hessian[i][j] += vals[k];
        if (i != j) {
          out.hessian[j][i] += vals[k];
        }
      }
    }
  }

  return out;
}

std::vector<double> eval_mvn_domain_grad(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma) {

  const GenzPack genz = pack_genz_ch(sigma, tape.perm);
  return eval_mvn_domain_grad_auto(tape, upper, mean, genz);
}

std::vector<double> eval_mvn_domain_grad_auto(
  MvnTape& tape,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const GenzPack& genz) {

  std::vector<double> analytic =
    analytic_mvn_domain_grad(upper, mean, genz, tape.lower);
  if (!analytic.empty()) {
    return analytic;
  }

  if (tape.value_only) {
    // Analytic path unavailable (singular / non-PD Sigma, bad scales, etc.).
    // Return NaNs so callers/optimizers can stop; do not throw under OpenMP.
    return std::vector<double>(
      tape.n_domain, std::numeric_limits<double>::quiet_NaN());
  }
  detail::init_qnorm_atomic();
  const CppAD::vector<double> x = detail::to_cppad_vector(pack_domain(upper, mean, genz));

  CppAD::sparse_rc<CppAD::vector<size_t>> jac_pattern;
  tape.fun.sparse_jac_rev(x, tape.pattern_grad, jac_pattern, detail::kJacColor, tape.work_grad);

  std::vector<double> grad_full(tape.n_domain, 0.0);
  const auto& cols = tape.pattern_grad.col();
  const auto& vals = tape.pattern_grad.val();
  for (size_t k = 0; k < tape.pattern_grad.nnz(); ++k) {
    grad_full[cols[k]] += vals[k];
  }
  return grad_full;
}

}  // namespace admvn
