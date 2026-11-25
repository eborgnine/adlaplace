#ifndef ADLAPLACE_ADFUN_HPP
#define ADLAPLACE_ADFUN_HPP

#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <mutex>
#include <string>
#include <vector>

#include "adlaplace/density_data.hpp"
#include "adlaplace/backend.hpp"
#include "adlaplace/parallel_build.hpp"
#include "adlaplace/rviews.hpp"

// Per-phase wall-clock breakdown of one adpack_sparsity() call (seconds),
// plus the colour counts the calibration sweeps report. Filled only when
// the caller passes a pointer; used for build profiling under verbose.
struct AdpackBuildTiming {
  double tape = 0.0;
  double discover_grad = 0.0;
  double discover_hess = 0.0;
  double build_patterns = 0.0;
  double jac = 0.0;
  double jac_inner = 0.0;
  double hes = 0.0;
  double hes_inner = 0.0;
  std::size_t n_color_jac = 0;
  std::size_t n_color_jac_inner = 0;
  std::size_t n_color_hes = 0;
  std::size_t n_color_hes_inner = 0;

  void add(const AdpackBuildTiming &o) {
    tape += o.tape;
    discover_grad += o.discover_grad;
    discover_hess += o.discover_hess;
    build_patterns += o.build_patterns;
    jac += o.jac;
    jac_inner += o.jac_inner;
    hes += o.hes;
    hes_inner += o.hes_inner;
    n_color_jac += o.n_color_jac;
    n_color_jac_inner += o.n_color_jac_inner;
    n_color_hes += o.n_color_hes;
    n_color_hes_inner += o.n_color_hes_inner;
  }
  double sparsity_total() const {
    return discover_grad + discover_hess + build_patterns + jac + jac_inner +
           hes + hes_inner;
  }
};

inline double adpack_timing_now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// Seed vector over the current tape domain (model.n_tape / tape_to_global).
inline CPPAD_TESTVECTOR(double)
    make_ad_params_seed(const Config &cfg, const density_data &model) {

  CPPAD_TESTVECTOR(double) ad_params_G(model.n_tape);
  for (std::size_t k = 0; k < model.n_tape; ++k) {
    const std::size_t g = model.tape_to_global[k];
    if (g < model.num_beta) {
      ad_params_G[k] = cfg.beta[g];
    } else if (g < model.num_beta + model.num_gamma) {
      ad_params_G[k] = cfg.gamma[g - model.num_beta];
    } else {
      ad_params_G[k] = cfg.theta[g - model.num_beta - model.num_gamma];
    }
  }
  return ad_params_G;
}

static const std::string JAC_COLOR = "cppad";
static const std::string HESS_COLOR = "cppad.symmetric";

template <class SizeVec, class ValVec>
inline ValVec &rcv_val(CppAD::sparse_rcv<SizeVec, ValVec> &rcv) {
  return const_cast<ValVec &>(rcv.val());
}

inline const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &empty_sparse_rc() {
  static CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> empty;
  return empty;
}

inline void
adpack_discover_grad(AdTape &gp,
                     CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &grad,
                     const std::size_t n_params) {

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> pattern_in;
  pattern_in.resize(n_params, n_params, n_params);
  for (std::size_t j = 0; j < n_params; ++j) {
    pattern_in.set(j, j, j);
  }

  const bool transpose = false;
  const bool dependency = false;
  const bool internal_bool = false;

  gp.fun.for_jac_sparsity(pattern_in, transpose, dependency, internal_bool,
                          grad);
}

inline void
adpack_discover_hessian(AdTape &gp,
                        CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian,
                        const std::size_t n_params) {

  CppAD::vectorBool select_domain(n_params), select_range(1);
  for (std::size_t j = 0; j < n_params; ++j) {
    select_domain[j] = true;
  }
  select_range[0] = true;

  const bool internal_bool = false;
  gp.fun.for_hes_sparsity(select_domain, select_range, internal_bool, hessian);
}

inline void adpack_attach_tape_maps(AdTape &gp, const density_data &model) {
  gp.n_global = model.num_full;
  gp.tape_to_global = model.tape_to_global;
  gp.gamma_row_to_tape.assign(model.num_gamma, model.n_tape);
  for (std::size_t k = 0; k < model.n_tape; ++k) {
    const std::size_t g = model.tape_to_global[k];
    if (g >= model.num_beta && g < model.num_beta + model.num_gamma) {
      gp.gamma_row_to_tape[g - model.num_beta] = k;
    }
  }
}

inline void adpack_fill_global_patterns(AdTape &gp) {
  const auto &map = gp.tape_to_global;
  const std::size_t n_tape = map.size();

  auto to_global = [&](std::size_t local) -> std::size_t {
    if (local >= n_tape) {
      Rcpp::stop("adpack_fill_global_patterns: local index out of range");
    }
    return map[local];
  };

  {
    const std::size_t nnz = gp.pattern_grad.nnz();
    const auto &cols = gp.pattern_grad.col();
    gp.grad_cols_global.resize(nnz);
    for (std::size_t d = 0; d < nnz; ++d) {
      gp.grad_cols_global[d] = to_global(cols[d]);
    }
  }
  {
    const std::size_t nnz = gp.pattern_grad_inner.nnz();
    const auto &cols = gp.pattern_grad_inner.col();
    gp.grad_inner_cols_global.resize(nnz);
    for (std::size_t d = 0; d < nnz; ++d) {
      gp.grad_inner_cols_global[d] = to_global(cols[d]);
    }
  }
  {
    const std::size_t nnz = gp.pattern_hessian.nnz();
    const auto &rows = gp.pattern_hessian.row();
    const auto &cols = gp.pattern_hessian.col();
    gp.hes_rows_global.resize(nnz);
    gp.hes_cols_global.resize(nnz);
    for (std::size_t d = 0; d < nnz; ++d) {
      gp.hes_rows_global[d] = to_global(rows[d]);
      gp.hes_cols_global[d] = to_global(cols[d]);
    }
  }
  {
    const std::size_t nnz = gp.pattern_hessian_inner.nnz();
    const auto &rows = gp.pattern_hessian_inner.row();
    const auto &cols = gp.pattern_hessian_inner.col();
    gp.hes_inner_rows_global.resize(nnz);
    gp.hes_inner_cols_global.resize(nnz);
    for (std::size_t d = 0; d < nnz; ++d) {
      gp.hes_inner_rows_global[d] = to_global(rows[d]);
      gp.hes_inner_cols_global[d] = to_global(cols[d]);
    }
  }
}

// Carrier for the patterns the calibration sweeps need. Hessian members keep
// both triangles; the upper-triangle views are what AdTape stores.
struct AdpackPatterns {
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad_inner;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_upper;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_inner_upper;
};

// Pure index work: inner/upper subsets of the given patterns. Touches no tape.
inline AdpackPatterns
adpack_build_patterns(AdTape &gp, const std::size_t n_params,
                      const std::vector<int> &subset,
                      const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &grad,
                      const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
                          &hessian_full,
                      const bool verbose = false) {

  gp.w.resize(1);
  gp.w[0] = 1.0;
  gp.x.resize(n_params);

  AdpackPatterns pats;
  pats.grad = grad;
  pats.hessian = hessian_full;

  std::vector<unsigned char> insubset(n_params, 0);
  for (std::size_t v : subset) {
    if (v < n_params) {
      insubset[v] = 1;
    }
  }

  const auto &full_cols = pats.grad.col();
  {
    std::size_t k = 0;
    for (std::size_t j : full_cols) {
      k += insubset[j];
    }
    pats.grad_inner.resize(1, n_params, k);

    std::size_t t = 0;
    for (std::size_t j : full_cols) {
      if (insubset[j])
        pats.grad_inner.set(t++, 0, j);
    }
  }

  {
    const auto &r = pats.hessian.row();
    const auto &c = pats.hessian.col();
    const std::size_t k = pats.hessian.nnz();

    std::size_t k_upper = 0, k_inner = 0, k_inner_upper = 0;

    for (std::size_t idx = 0; idx < k; ++idx) {
      const std::size_t i = r[idx], j = c[idx];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper)
        ++k_upper;
      if (is_inner)
        ++k_inner;
      if (is_upper && is_inner)
        ++k_inner_upper;
    }

    pats.hessian_upper.resize(n_params, n_params, k_upper);
    pats.hessian_inner.resize(n_params, n_params, k_inner);
    pats.hessian_inner_upper.resize(n_params, n_params, k_inner_upper);

    std::size_t tu = 0, ti = 0, tiu = 0;
    for (std::size_t idx = 0; idx < k; ++idx) {
      const std::size_t i = r[idx], j = c[idx];
      const bool is_upper = (i <= j);
      const bool is_inner = insubset[i] && insubset[j];

      if (is_upper)
        pats.hessian_upper.set(tu++, i, j);
      if (is_inner)
        pats.hessian_inner.set(ti++, i, j);
      if (is_upper && is_inner)
        pats.hessian_inner_upper.set(tiu++, i, j);
    }
  }

  if (verbose) {
    Rcpp::Rcout << "  sparsity: grad " << pats.grad.nnz() << ", grad_inner "
                << pats.grad_inner.nnz() << ", hes full" << pats.hessian.nnz()
                << ", hes upper" << pats.hessian_upper.nnz() << ", hes_inner "
                << pats.hessian_inner_upper.nnz() << "\n";
  }

  return pats;
}

// Assign the four rcv patterns onto the pack and map them to global indices.
// Safe with or without a recorded tape.
inline void adpack_install_patterns(AdTape &gp, const AdpackPatterns &pats) {
  gp.pattern_grad =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.grad);
  gp.pattern_grad_inner =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.grad_inner);
  gp.pattern_hessian =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.hessian_upper);
  gp.pattern_hessian_inner =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.hessian_inner_upper);
  adpack_fill_global_patterns(gp);
}

// One-time CppAD coloring / reverse sweeps that fill the work_* objects
// EvalShard later reuses. Requires a recorded tape. Interleaves rcv
// construction with each sweep, matching the original adpack_sparsity order.
inline void adpack_calibrate(const CPPAD_TESTVECTOR(double) & x, AdTape &gp,
                             AdpackPatterns &pats,
                             AdpackBuildTiming *timing = nullptr) {
  gp.fun.Forward(0, x);

  gp.pattern_grad =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.grad);
  {
    const double t0 = timing ? adpack_timing_now() : 0.0;
    const std::size_t nc = gp.fun.sparse_jac_rev(
        x, gp.pattern_grad, pats.grad, JAC_COLOR, gp.work_grad);
    if (timing) {
      timing->jac += adpack_timing_now() - t0;
      timing->n_color_jac += nc;
    }
  }

  gp.pattern_grad_inner =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.grad_inner);
  {
    const double t0 = timing ? adpack_timing_now() : 0.0;
    const std::size_t nc = gp.fun.sparse_jac_rev(
        x, gp.pattern_grad_inner, pats.grad_inner, JAC_COLOR,
        gp.work_inner_grad);
    if (timing) {
      timing->jac_inner += adpack_timing_now() - t0;
      timing->n_color_jac_inner += nc;
    }
  }

  gp.pattern_hessian =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.hessian_upper);
  gp.pattern_hessian_inner =
      CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)>(
          pats.hessian_inner_upper);

  // Grad-only packs leave Hessian patterns empty; skip sparse_hes calibration.
  if (pats.hessian.nnz() > 0 || pats.hessian_upper.nnz() > 0) {
    const double t0 = timing ? adpack_timing_now() : 0.0;
    const std::size_t nc = gp.fun.sparse_hes(
        x, gp.w, gp.pattern_hessian, pats.hessian, HESS_COLOR, gp.work_hess);
    if (timing) {
      timing->hes += adpack_timing_now() - t0;
      timing->n_color_hes += nc;
    }
    const double t1 = timing ? adpack_timing_now() : 0.0;
    const std::size_t nc_inner = gp.fun.sparse_hes(
        x, gp.w, gp.pattern_hessian_inner, pats.hessian_inner, HESS_COLOR,
        gp.work_inner_hess);
    if (timing) {
      timing->hes_inner += adpack_timing_now() - t1;
      timing->n_color_hes_inner += nc_inner;
    }
  }

  adpack_fill_global_patterns(gp);
  // Drop the calibration Taylor buffer in the same DSO that recorded the
  // tape. Later drain from adlaplace.so then sees cap_order 0 and is a no-op.
  gp.fun.capacity_order(0);
}

inline void adpack_sparsity(const CPPAD_TESTVECTOR(double) & x,
                            const std::vector<int> &subset, AdTape &gp,
                            const bool verbose = false,
                            const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
                                &hessian = empty_sparse_rc(),
                            const bool hessian_sparsity = true,
                            AdpackBuildTiming *timing = nullptr) {

  const std::size_t n_params = x.size();

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> grad;
  {
    const double t0 = timing ? adpack_timing_now() : 0.0;
    adpack_discover_grad(gp, grad, n_params);
    if (timing) {
      timing->discover_grad += adpack_timing_now() - t0;
    }
  }

  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hessian_here;
  if (hessian.nnz() > 0) {
    if (verbose) {
      Rcpp::Rcout << "using external Hessian sparsity pattern\n";
    }
    hessian_here = hessian;
  } else if (hessian_sparsity) {
    if (verbose) {
      Rcpp::Rcout << "discovering Hessian sparsity pattern\n";
    }
    const double t0 = timing ? adpack_timing_now() : 0.0;
    adpack_discover_hessian(gp, hessian_here, n_params);
    if (timing) {
      timing->discover_hess += adpack_timing_now() - t0;
    }
  } else {
    if (verbose) {
      Rcpp::Rcout << "skipping Hessian sparsity (grad-only)\n";
    }
    hessian_here = empty_sparse_rc();
  }

  AdpackPatterns pats;
  {
    const double t0 = timing ? adpack_timing_now() : 0.0;
    pats = adpack_build_patterns(gp, n_params, subset, grad, hessian_here,
                                 verbose);
    if (timing) {
      timing->build_patterns += adpack_timing_now() - t0;
    }
  }
  adpack_calibrate(x, gp, pats, timing);
}

inline size_t count_obs_shards(const density_data &model, const Rcpp::List &config) {
  const Config cfg(config);
  size_t ng = cfg.obs_groups.ncol();
  if (ng == 0) {
    if (model.y.size() > 0) {
      ng = 1;
    } else if (model.elgm_matrix.p.size() > 0) {
      ng = static_cast<size_t>(model.elgm_matrix.ncol());
    }
  }
  return ng > 0 ? ng : 1;
}

// Defined in obs_sparsity.hpp (included below). Forward-declared so
// build_ad_fun_obs compiles if that header's include is skipped by
// include-order (obs_sparsity.hpp -> ad_pack_random.hpp -> ad_pack.hpp).
inline CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)>
obs_hessian_sparsity(const density_data &model, const Config &cfg,
                     std::size_t Dgroup);

#include "adlaplace/obs_sparsity.hpp"

inline std::vector<AdTape> build_ad_fun_obs(const density_data &model_in,
                                               const Rcpp::List &config,
                                               LogDensObsFn log_dens) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model_in);
  const size_t ng = count_obs_shards(model_in, config);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_obs groups " << ng << "\n";
  }

  std::vector<AdTape> result(ng);
  std::vector<AdpackBuildTiming> timings(ng);
  AdpackBuildTiming *timing_ptr = cfg.verbose ? timings.data() : nullptr;

  const std::size_t n_build_threads =
      adlaplace_parallel_build::effective_threads(cfg.num_threads, ng);

  if (n_build_threads <= 1) {
    // Serial path: unchanged two-pass build (tape all groups, then
    // sparsity/calibration for all groups).
    density_data model = model_in;

    for (size_t d = 0; d < ng; ++d) {
      if (cfg.verbose) {
        Rcpp::Rcout << "  taping observation group " << (d + 1) << " / " << ng
                    << "\n";
      }
      model.apply_tape_domain(cfg, "obs", d);
      const CPPAD_TESTVECTOR(double) ad_params_G =
          make_ad_params_seed(cfg, model);
      const double t0 = timing_ptr ? adpack_timing_now() : 0.0;
      CppAD::vector<CppAD::AD<double>> ad_params(model.n_tape);
      for (size_t j = 0; j < model.n_tape; ++j) {
        ad_params[j] = ad_params_G[j];
      }
      CppAD::Independent(ad_params);
      auto result_here = log_dens(ad_params, model, cfg, d);
      CppAD::ADFun<double> fun(ad_params, result_here);
      if (timing_ptr) {
        timings[d].tape += adpack_timing_now() - t0;
      }
      result[d].fun = std::move(fun);
      adpack_attach_tape_maps(result[d], model);
    }

    for (size_t d = 0; d < ng; ++d) {
      if (cfg.verbose) {
        Rcpp::Rcout << "  sparsity observation group " << (d + 1) << " / "
                    << ng << "\n";
      }
      model.apply_tape_domain(cfg, "obs", d);
      const CPPAD_TESTVECTOR(double) ad_params_G =
          make_ad_params_seed(cfg, model);
      CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hess = empty_sparse_rc();
      if (cfg.hessian_sparsity && use_analytic_obs_hessian(cfg)) {
        hess = obs_hessian_sparsity(model, cfg, d);
      }
      adpack_sparsity(ad_params_G, model.seq_gamma, result[d], cfg.verbose,
                      hess, cfg.hessian_sparsity,
                      timing_ptr ? &timings[d] : nullptr);
    }
  } else {
    // Parallel path: groups are independent (own ADFun, patterns and
    // sparse work objects), so tape + calibrate one group per task.
    // Shard costs are heterogeneous by construction, hence dynamic
    // scheduling. No R API calls inside the region: verbose is silenced
    // for workers and any error is recorded and rethrown serially.
    Config cfg_worker = cfg;
    cfg_worker.verbose = false;

    // Per-group tape domains, prepared serially: density_data copies
    // hold Rcpp vectors, and copying/destroying those on OpenMP workers
    // corrupts R's (main-thread) protect stack. Workers therefore only
    // read models[d]; all Rcpp object lifetimes stay on the main thread.
    std::vector<density_data> models;
    models.reserve(ng);
    for (std::size_t d = 0; d < ng; ++d) {
      models.push_back(model_in);
      models.back().apply_tape_domain(cfg, "obs", d);
    }

    std::atomic<bool> failed(false);
    std::mutex error_mutex;
    std::string first_error;

    {
      adlaplace_parallel_build::Guard parallel_guard(n_build_threads);

#pragma omp parallel for schedule(dynamic, 1) num_threads(static_cast<int>(n_build_threads))
      for (std::ptrdiff_t di = 0; di < static_cast<std::ptrdiff_t>(ng);
           ++di) {
        if (failed.load(std::memory_order_relaxed)) {
          continue;
        }
        const std::size_t d = static_cast<std::size_t>(di);
        try {
          const density_data &model = models[d];
          const CPPAD_TESTVECTOR(double) ad_params_G =
              make_ad_params_seed(cfg, model);
          AdpackBuildTiming *tp = timing_ptr ? &timings[d] : nullptr;
          const double t0 = tp ? adpack_timing_now() : 0.0;
          CppAD::vector<CppAD::AD<double>> ad_params(model.n_tape);
          for (size_t j = 0; j < model.n_tape; ++j) {
            ad_params[j] = ad_params_G[j];
          }
          CppAD::Independent(ad_params);
          auto result_here = log_dens(ad_params, model, cfg_worker, d);
          CppAD::ADFun<double> fun(ad_params, result_here);
          if (tp) {
            tp->tape += adpack_timing_now() - t0;
          }
          result[d].fun = std::move(fun);
          adpack_attach_tape_maps(result[d], model);

          CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> hess =
              empty_sparse_rc();
          if (cfg.hessian_sparsity && use_analytic_obs_hessian(cfg)) {
            hess = obs_hessian_sparsity(model, cfg, d);
          }
          adpack_sparsity(ad_params_G, model.seq_gamma, result[d], false,
                          hess, cfg.hessian_sparsity, tp);
        } catch (const std::exception &e) {
          // If the failure interrupted taping, abandon this thread's
          // half-recorded tape; a no-op when not recording.
          CppAD::AD<double>::abort_recording();
          std::lock_guard<std::mutex> lock(error_mutex);
          if (first_error.empty()) {
            first_error = e.what();
          }
          failed.store(true, std::memory_order_relaxed);
        } catch (...) {
          CppAD::AD<double>::abort_recording();
          std::lock_guard<std::mutex> lock(error_mutex);
          if (first_error.empty()) {
            first_error = "unknown error during parallel observation build";
          }
          failed.store(true, std::memory_order_relaxed);
        }
      }
    } // parallel_guard teardown: restore sequential CppAD mode

    if (failed.load(std::memory_order_relaxed)) {
      Rcpp::stop("build_ad_fun_obs (parallel): %s", first_error.c_str());
    }
  }

  if (cfg.verbose) {
    AdpackBuildTiming total;
    for (const auto &t : timings) {
      total.add(t);
    }
    Rcpp::Rcout << "build_ad_fun_obs timing: threads " << n_build_threads
                << ", tape " << total.tape << "s, discover_grad "
                << total.discover_grad << "s, discover_hess "
                << total.discover_hess << "s, build_patterns "
                << total.build_patterns << "s, jac " << total.jac
                << "s (colors " << total.n_color_jac << "), jac_inner "
                << total.jac_inner << "s (colors "
                << total.n_color_jac_inner << "), hes " << total.hes
                << "s (colors " << total.n_color_hes << "), hes_inner "
                << total.hes_inner << "s (colors "
                << total.n_color_hes_inner << ")\n";
  }

  return result;
}

inline AdTape build_ad_fun_parameters(
    const density_data &model_in, const Rcpp::List &config, LogDensSingleDataFn log_dens,
    const CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> &hessian = empty_sparse_rc()) {

  const Config cfg(config);
  validate_config_matches_model(cfg, model_in);

  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_parameters: taping...\n";
  }

  density_data model = model_in;
  model.apply_tape_domain(cfg, "all", 0);

  const CPPAD_TESTVECTOR(double) ad_params_G = make_ad_params_seed(cfg, model);

  CppAD::vector<CppAD::AD<double>> ad_params(model.n_tape);
  for (size_t d = 0; d < model.n_tape; ++d) {
    ad_params[d] = ad_params_G[d];
  }

  CppAD::Independent(ad_params);

  CppAD::vector<CppAD::AD<double>> result_here =
      log_dens(ad_params, model, cfg);

  CppAD::ADFun<double> fun(ad_params, result_here);

  AdTape pack;
  pack.fun = std::move(fun);
  pack.owner_thread_assigned = false;
  adpack_attach_tape_maps(pack, model);
  if (cfg.verbose) {
    Rcpp::Rcout << "build_ad_fun_parameters: computing sparsity...\n";
  }
  adpack_sparsity(ad_params_G, model.seq_gamma, pack, cfg.verbose, hessian,
                  cfg.hessian_sparsity);
  return pack;
}

#endif
