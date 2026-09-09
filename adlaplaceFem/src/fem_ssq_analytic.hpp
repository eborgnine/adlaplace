#ifndef ADLAPLACEFEM_FEM_SSQ_ANALYTIC_HPP
#define ADLAPLACEFEM_FEM_SSQ_ANALYTIC_HPP

// Closed-form derivatives of the FEM quadratic form
//
//   f(gamma, theta) = -0.5 * gamma' Q(theta) gamma,
//   Q(theta)        = sum_j c_j(theta) * M_j,   M = (C, G, G2, [G3]).
//
// Q is a linear combination of constant Grams whose m = 3 or 4 coefficients
// depend only on the two thetas, and f is quadratic in gamma, so every
// derivative the Laplace approximation needs is available in closed form:
//
//   df/dgamma            = -Q gamma
//   df/dtheta_k          = -0.5 sum_j c'_jk  (gamma' M_j gamma)
//   d2f/dgamma dgamma    = -Q                (the combined values themselves)
//   d2f/dgamma dtheta_k  = -(dQ/dtheta_k) gamma
//   d2f/dtheta_j dtheta_k= -0.5 sum_l c''_ljk (gamma' M_l gamma)
//
// Taping the quadratic form instead puts every nonzero of Q on the tape, and
// sparse_hes then sweeps it once per color while trace_hinv_t sweeps it four
// times per assigned column. All of that collapses to a handful of O(nnz)
// passes here.
//
// The third-order trace direction is scattered into gamma slots only, so the
// order-2 Taylor coefficient is -0.5 * v' Q(theta) v and its gradient is
// identically zero on every gamma coordinate, leaving only
//   d/dtheta_k = -0.5 * v' (dQ/dtheta_k) v.

#include <cppad/cppad.hpp>

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace femssq {

using PatternRcv =
    CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>;

// Number of thetas this density uses (range, sd).
constexpr std::size_t kNumTheta = 2;

// Immutable per-model data, shared by every clone of the shard.
struct Data {
  int alpha = 2;
  // Rows of Q, which is also the number of gamma.
  std::size_t n = 0;
  // Upper-triangle CSC pattern of Q, from fem_precision_payload().
  std::vector<int> Q_p;
  std::vector<int> Q_i;
  // 1 on the diagonal, 2 off it, so one pass gives a full quadratic form.
  std::vector<double> w;
  // Gram value vectors aligned with Q_p/Q_i, one per coefficient.
  std::vector<std::vector<double>> M;
  // Gamma row -> tape index.
  std::vector<std::size_t> gidx;
  // Tape indices of the two thetas.
  std::size_t t_tape[kNumTheta] = {0, 0};
  // Domain 2, Range m map theta -> coefficients. Never evaluated through this
  // handle; each shard takes its own copy so threads never share a tape.
  CppAD::ADFun<double> coef_fun;

  std::size_t m() const { return M.size(); }
  std::size_t nnz() const { return Q_i.size(); }
};

// ---------------------------------------------------------------------------
// Primitives over the upper-triangle CSC storage
// ---------------------------------------------------------------------------

// out[k] = sum_j coef[j] * M_j[k], the values of a single combined matrix.
inline void combine(const Data &d, const double *coef,
                    std::vector<double> &out) {
  const std::size_t nnz = d.nnz();
  out.assign(nnz, 0.0);
  for (std::size_t j = 0; j < d.m(); ++j) {
    const double cj = coef[j];
    if (cj == 0.0) {
      continue;
    }
    const std::vector<double> &Mj = d.M[j];
    for (std::size_t k = 0; k < nnz; ++k) {
      out[k] += cj * Mj[k];
    }
  }
}

// y = A * x for a symmetric A stored as upper-triangle CSC.
inline void sym_matvec(const Data &d, const std::vector<double> &A,
                       const std::vector<double> &x, std::vector<double> &y) {
  y.assign(d.n, 0.0);
  for (std::size_t col = 0; col < d.n; ++col) {
    const double xc = x[col];
    // Mirrored contributions land on y[col]; accumulate them separately so
    // the diagonal entry is not counted twice.
    double mirrored = 0.0;
    for (int pos = d.Q_p[col]; pos < d.Q_p[col + 1]; ++pos) {
      const std::size_t k = static_cast<std::size_t>(pos);
      const std::size_t row = static_cast<std::size_t>(d.Q_i[k]);
      y[row] += A[k] * xc;
      if (row != col) {
        mirrored += A[k] * x[row];
      }
    }
    y[col] += mirrored;
  }
}

// out[j] = x' M_j x for every Gram, in one pass over the shared pattern.
inline void sym_quad_all(const Data &d, const std::vector<double> &x,
                         std::vector<double> &out) {
  const std::size_t m = d.m();
  out.assign(m, 0.0);
  for (std::size_t col = 0; col < d.n; ++col) {
    const double xc = x[col];
    if (xc == 0.0) {
      continue;
    }
    for (int pos = d.Q_p[col]; pos < d.Q_p[col + 1]; ++pos) {
      const std::size_t k = static_cast<std::size_t>(pos);
      const double fac = d.w[k] * xc * x[static_cast<std::size_t>(d.Q_i[k])];
      if (fac == 0.0) {
        continue;
      }
      for (std::size_t j = 0; j < m; ++j) {
        out[j] += fac * d.M[j][k];
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Pattern slot maps
// ---------------------------------------------------------------------------

// Tape index -> gamma row, or -1.
inline std::vector<int> tape_to_gamma_row(const Data &d, std::size_t n_tape) {
  std::vector<int> out(n_tape, -1);
  for (std::size_t r = 0; r < d.gidx.size(); ++r) {
    if (d.gidx[r] < n_tape) {
      out[d.gidx[r]] = static_cast<int>(r);
    }
  }
  return out;
}

// (row, col) of the upper triangle -> position in Q_i.
inline std::unordered_map<std::uint64_t, int> build_q_index(const Data &d) {
  std::unordered_map<std::uint64_t, int> map;
  map.reserve(d.nnz() * 2);
  for (std::size_t col = 0; col < d.n; ++col) {
    for (int pos = d.Q_p[col]; pos < d.Q_p[col + 1]; ++pos) {
      const std::uint64_t key =
          static_cast<std::uint64_t>(d.Q_i[static_cast<std::size_t>(pos)]) +
          static_cast<std::uint64_t>(col) * static_cast<std::uint64_t>(d.n);
      map.emplace(key, pos);
    }
  }
  return map;
}

// What each nonzero of a gradient pattern is: >= 0 is a gamma row, -1 - k is
// theta k, and kIgnore is a slot with no analytic contribution.
struct GradSlots {
  static constexpr int kIgnore = -100;
  std::vector<int> target;
};

inline GradSlots build_grad_slots(const Data &d, const PatternRcv &pat,
                                  const std::vector<int> &t2g) {
  GradSlots s;
  const std::size_t nnz = pat.nnz();
  s.target.assign(nnz, GradSlots::kIgnore);
  const auto &cols = pat.col();
  for (std::size_t k = 0; k < nnz; ++k) {
    const std::size_t j = cols[k];
    bool done = false;
    for (std::size_t t = 0; t < kNumTheta && !done; ++t) {
      if (j == d.t_tape[t]) {
        s.target[k] = -1 - static_cast<int>(t);
        done = true;
      }
    }
    if (!done && j < t2g.size() && t2g[j] >= 0) {
      s.target[k] = t2g[j];
    }
  }
  return s;
}

// What each nonzero of a Hessian pattern is.
struct HesSlots {
  enum Kind : unsigned char {
    kGammaGamma = 0,  // a = position in Q_i
    kGammaTheta = 1,  // a = gamma row, b = theta index
    kThetaTheta = 2,  // a, b = theta indices
    kIgnore = 3
  };
  std::vector<unsigned char> kind;
  std::vector<int> a;
  std::vector<int> b;
};

inline HesSlots
build_hes_slots(const Data &d, const PatternRcv &pat,
                const std::vector<int> &t2g,
                const std::unordered_map<std::uint64_t, int> &qmap) {
  HesSlots s;
  const std::size_t nnz = pat.nnz();
  s.kind.assign(nnz, HesSlots::kIgnore);
  s.a.assign(nnz, 0);
  s.b.assign(nnz, 0);
  const auto &rows = pat.row();
  const auto &cols = pat.col();

  const auto theta_of = [&d](std::size_t tape) -> int {
    for (std::size_t t = 0; t < kNumTheta; ++t) {
      if (tape == d.t_tape[t]) {
        return static_cast<int>(t);
      }
    }
    return -1;
  };

  for (std::size_t k = 0; k < nnz; ++k) {
    const std::size_t i = rows[k];
    const std::size_t j = cols[k];
    const int ri = (i < t2g.size()) ? t2g[i] : -1;
    const int rj = (j < t2g.size()) ? t2g[j] : -1;
    const int ti = theta_of(i);
    const int tj = theta_of(j);

    if (ri >= 0 && rj >= 0) {
      const std::uint64_t lo = static_cast<std::uint64_t>(ri < rj ? ri : rj);
      const std::uint64_t hi = static_cast<std::uint64_t>(ri < rj ? rj : ri);
      const auto it =
          qmap.find(lo + hi * static_cast<std::uint64_t>(d.n));
      if (it != qmap.end()) {
        s.kind[k] = HesSlots::kGammaGamma;
        s.a[k] = it->second;
      }
    } else if (ri >= 0 && tj >= 0) {
      s.kind[k] = HesSlots::kGammaTheta;
      s.a[k] = ri;
      s.b[k] = tj;
    } else if (rj >= 0 && ti >= 0) {
      s.kind[k] = HesSlots::kGammaTheta;
      s.a[k] = rj;
      s.b[k] = ti;
    } else if (ti >= 0 && tj >= 0) {
      s.kind[k] = HesSlots::kThetaTheta;
      s.a[k] = ti;
      s.b[k] = tj;
    }
  }
  return s;
}

// ---------------------------------------------------------------------------
// Per-instance evaluation state
// ---------------------------------------------------------------------------

struct Workspace {
  // Own copy of the coefficient tape: CppAD sweeps are not thread safe, and
  // clones run on different threads.
  CppAD::ADFun<double> coef_fun;

  // Coefficient values and derivatives at the current theta.
  std::vector<double> c;                 // m
  std::vector<double> dc;                // m * kNumTheta, dc[j * nt + k]
  std::vector<double> d2c;               // m * nt * nt
  std::vector<double> dc_col;            // scratch column of dc

  std::vector<double> Q_x;               // nnz, sum_j c_j M_j
  std::vector<double> dQ_x[kNumTheta];   // nnz, sum_j c'_jk M_j

  std::vector<double> gamma;             // n
  std::vector<double> Qg;                // n
  std::vector<double> dQg[kNumTheta];    // n
  std::vector<double> quad;              // m, gamma' M_j gamma

  // trace_hinv_t scratch: a stamped dense direction avoids clearing n
  // entries per column.
  std::vector<double> v_dense;
  std::vector<int> stamp;
  std::vector<int> touched;
  int stamp_tag = 0;

  void init(const Data &d) {
    coef_fun = d.coef_fun;
    v_dense.assign(d.n, 0.0);
    stamp.assign(d.n, 0);
    stamp_tag = 0;
  }
};

// Evaluate the coefficient map and build the combined value vectors.
// `order` 0 gives the value only, 1 adds dc and dQ/dtheta, 2 adds d2c. The
// coefficient tape has Domain 2, so its sweeps cost nothing next to the
// O(m * nnz) combine passes.
inline void update_theta(const Data &d, Workspace &ws, const double *theta,
                         int order) {
  const std::size_t m = d.m();
  const std::size_t nt = kNumTheta;

  CPPAD_TESTVECTOR(double) tv(nt);
  for (std::size_t k = 0; k < nt; ++k) {
    tv[k] = theta[k];
  }

  const CPPAD_TESTVECTOR(double) cv = ws.coef_fun.Forward(0, tv);
  ws.c.assign(m, 0.0);
  for (std::size_t j = 0; j < m; ++j) {
    ws.c[j] = cv[j];
  }
  combine(d, ws.c.data(), ws.Q_x);
  if (order < 1) {
    return;
  }

  const CPPAD_TESTVECTOR(double) jac = ws.coef_fun.Jacobian(tv);
  ws.dc.assign(m * nt, 0.0);
  for (std::size_t j = 0; j < m; ++j) {
    for (std::size_t k = 0; k < nt; ++k) {
      ws.dc[j * nt + k] = jac[j * nt + k];
    }
  }

  ws.dc_col.assign(m, 0.0);
  for (std::size_t k = 0; k < nt; ++k) {
    for (std::size_t j = 0; j < m; ++j) {
      ws.dc_col[j] = ws.dc[j * nt + k];
    }
    combine(d, ws.dc_col.data(), ws.dQ_x[k]);
  }
  if (order < 2) {
    return;
  }

  ws.d2c.assign(m * nt * nt, 0.0);
  for (std::size_t j = 0; j < m; ++j) {
    const CPPAD_TESTVECTOR(double) hes = ws.coef_fun.Hessian(tv, j);
    for (std::size_t p = 0; p < nt * nt; ++p) {
      ws.d2c[j * nt * nt + p] = hes[p];
    }
  }
}

// -0.5 * gamma' Q gamma, given the quadratic forms already computed.
inline double value_from_quad(const Data &d, const Workspace &ws) {
  double acc = 0.0;
  for (std::size_t j = 0; j < d.m(); ++j) {
    acc += ws.c[j] * ws.quad[j];
  }
  return -0.5 * acc;
}

// df/dtheta_k = -0.5 sum_j c'_jk (gamma' M_j gamma)
inline double grad_theta(const Data &d, const Workspace &ws, std::size_t k) {
  double acc = 0.0;
  for (std::size_t j = 0; j < d.m(); ++j) {
    acc += ws.dc[j * kNumTheta + k] * ws.quad[j];
  }
  return -0.5 * acc;
}

// d2f/dtheta_a dtheta_b = -0.5 sum_j c''_jab (gamma' M_j gamma)
inline double hes_theta(const Data &d, const Workspace &ws, std::size_t a,
                        std::size_t b) {
  double acc = 0.0;
  for (std::size_t j = 0; j < d.m(); ++j) {
    acc += ws.d2c[j * kNumTheta * kNumTheta + a * kNumTheta + b] * ws.quad[j];
  }
  return -0.5 * acc;
}

} // namespace femssq

#endif
