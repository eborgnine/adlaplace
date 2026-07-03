#include "mvn_tape.hpp"

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

template <typename ADType>
ADType qmc_integrand(
  const std::vector<ADType>& bs,
  const std::vector<double>& as,
  const std::vector<std::vector<double>>& ch,
  const std::vector<double>& w) {

  const std::size_t n = bs.size();
  const double ct0 = ch[0][0];
  const double ai0 = as[0];
  const ADType bi0 = bs[0];

  ADType c1 = ADType(0.0);
  if (ai0 > -kNine * ct0) {
    if (ai0 < kNine * ct0) {
      c1 = detail::pnorm_ad(ADType(ai0) / ct0);
    } else {
      c1 = ADType(1.0);
    }
  }

  ADType d1 = ADType(0.0);
  if (bi0 > ADType(-kNine * ct0)) {
    if (bi0 < ADType(kNine * ct0)) {
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
      s += ADType(ch[i][k]) * y[k];
    }

    const double ct = ch[i][i];
    const double asi = as[i];
    const ADType aicnt = ADType(asi) - s;
    const ADType bicnt = bs[i] - s;

    ADType ci = ADType(1.0);
    ADType di = ADType(1.0);

    if (aicnt < ADType(-kNine * ct)) {
      ci = ADType(0.0);
    } else if (CppAD::abs(aicnt) < ADType(kNine * ct)) {
      ci = detail::pnorm_ad(aicnt / ct);
    }

    if (bicnt < ADType(-kNine * ct)) {
      di = ADType(0.0);
    } else if (CppAD::abs(bicnt) < ADType(kNine * ct)) {
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
  const MvnSetup& setup) {

  const std::size_t n = setup.n;
  std::vector<double> bs(n);
  for (std::size_t i = 0; i < n; ++i) {
    const int j = setup.perm[i];
    bs[i] = (upper[j] - mean[j]) / setup.scale[j];
  }
  return bs;
}

double qmc_error(
  const MvnSetup& setup,
  const std::vector<double>& upper,
  const std::vector<double>& mean) {

  const std::vector<double> bs = make_bs(upper, mean, setup);
  double p = 0.0;
  double e = 0.0;

  for (std::size_t shift = 0; shift < setup.n_shifts; ++shift) {
    double sm = 0.0;
    for (std::size_t pt = 0; pt < setup.n_points; ++pt) {
      sm += qmc_integrand_double(bs, setup.as, setup.ch, setup.qmc_w[shift][pt]);
    }
    sm /= static_cast<double>(setup.n_points);
    const double dm = (sm - p) / static_cast<double>(shift + 1);
    p += dm;
    e = (static_cast<double>(shift) * e / static_cast<double>(shift + 1)) + dm * dm;
  }
  (void)p;
  return 3.0 * std::sqrt(e);
}

}  // namespace

MvnSetup prepare_mvn_setup(
  const std::vector<double>& lower,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  const std::vector<std::vector<double>>& sigma,
  std::size_t n_points,
  std::size_t n_shifts,
  unsigned int seed) {

  const std::size_t n = sigma.size();
  if (n == 0) {
    throw std::runtime_error("sigma must be at least 1x1");
  }

  MvnSetup setup;
  setup.n = n;
  setup.n_points = n_points;
  setup.n_shifts = n_shifts;
  setup.mean = mean;
  setup.perm.resize(n);
  setup.inv_perm.resize(n);
  setup.scale.resize(n);
  std::iota(setup.perm.begin(), setup.perm.end(), 0);

  std::vector<std::vector<double>> c = sigma;
  std::vector<double> ap(n);
  std::vector<double> bp(n);
  for (std::size_t i = 0; i < n; ++i) {
    setup.scale[i] = std::sqrt(sigma[i][i]);
    ap[i] = is_neg_inf(lower[i]) ? -std::numeric_limits<double>::infinity()
                                 : (lower[i] - mean[i]) / setup.scale[i];
    bp[i] = is_pos_inf(upper[i]) ? std::numeric_limits<double>::infinity()
                                 : (upper[i] - mean[i]) / setup.scale[i];
  }

  chlrdr(c, ap, bp, setup.perm);
  setup.as = ap;
  setup.ch = c;

  for (std::size_t i = 0; i < n; ++i) {
    setup.inv_perm[setup.perm[i]] = static_cast<int>(i);
  }

  const std::vector<double> q = primes_sqrt(n);
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> unif(0.0, 1.0);

  setup.qmc_w.assign(
    n_shifts,
    std::vector<std::vector<double>>(
      n_points,
      std::vector<double>(n > 1 ? n - 1 : 0, 0.0)));

  for (std::size_t shift = 0; shift < n_shifts; ++shift) {
    const double xr = unif(rng);
    for (std::size_t pt = 0; pt < n_points; ++pt) {
      const double cnt = static_cast<double>(pt + 1);
      for (std::size_t dim = 0; dim + 1 < n; ++dim) {
        double x = std::fmod(cnt * q[dim] + xr, 1.0);
        if (x < 0.0) {
          x += 1.0;
        }
        setup.qmc_w[shift][pt][dim] = std::abs(2.0 * x - 1.0);
      }
    }
  }

  return setup;
}

CppAD::ADFun<double> build_mvn_tape(const MvnSetup& setup) {
  detail::init_qnorm_atomic();

  const std::size_t n = setup.n;
  std::vector<CppAD::AD<double>> upper_ad(n);
  CppAD::Independent(upper_ad);

  std::vector<CppAD::AD<double>> bs(n);
  for (std::size_t i = 0; i < n; ++i) {
    const int j = setup.perm[i];
    bs[i] = (upper_ad[j] - setup.mean[j]) / setup.scale[j];
  }

  CppAD::AD<double> total = CppAD::AD<double>(0.0);
  const double norm = static_cast<double>(setup.n_shifts * setup.n_points);

  for (std::size_t shift = 0; shift < setup.n_shifts; ++shift) {
    for (std::size_t pt = 0; pt < setup.n_points; ++pt) {
      total += qmc_integrand(bs, setup.as, setup.ch, setup.qmc_w[shift][pt]);
    }
  }
  total /= norm;

  std::vector<CppAD::AD<double>> y(1);
  y[0] = total;
  CppAD::ADFun<double> fun(upper_ad, y);
  return fun;
}

MvnResult eval_mvn_tape(
  const MvnSetup& setup,
  const std::vector<double>& upper,
  const std::vector<double>& mean,
  CppAD::ADFun<double>* pre_taped) {

  const std::size_t n = setup.n;
  MvnResult out;
  out.gradient.assign(n, 0.0);
  out.hessian.assign(n, std::vector<double>(n, 0.0));

  detail::init_qnorm_atomic();

  CppAD::ADFun<double> owned;
  CppAD::ADFun<double>* fun = pre_taped;
  if (!fun) {
    owned = build_mvn_tape(setup);
    fun = &owned;
  }

  const std::vector<double> x = upper;
  out.value = (*fun).Forward(0, x)[0];
  out.gradient = (*fun).Jacobian(x);
  const std::vector<double> h = (*fun).Hessian(x, 0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      out.hessian[i][j] = h[i * n + j];
    }
  }
  out.error = qmc_error(setup, upper, mean);

  return out;
}

}  // namespace admvn
