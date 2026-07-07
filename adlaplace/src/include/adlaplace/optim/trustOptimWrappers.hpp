#ifndef ADLAPLACE_AD_FUNC_OPT_HPP
#define ADLAPLACE_AD_FUNC_OPT_HPP

#include <Rcpp.h>
#include <Rinternals.h>

#include <Eigen/Sparse>
#include <numeric>
#include <vector>

#include <omp.h>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/api/backend.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"
#include "adlaplace/runtime/thread_affinity_debug.hpp"

struct AD_Func_Opt {
  ad_fun *const shards;
  const std::vector<std::vector<size_t>> *const thread_groups;
  const bool inner;
  const int num_threads;

  const size_t Nparams;
  const size_t n_shards_backend;
  const size_t Nbeta_backend;
  const size_t Ngamma_backend;
  const size_t Ntheta_backend;
  const size_t nvars_opt;
  const size_t var_offset;

  const Eigen::SparseMatrix<int, Eigen::ColMajor, int> Htemplate;
  const std::vector<std::vector<int>> hess_maps;
  const std::vector<size_t> all_shards;

  std::vector<double> parameters;
  std::vector<double> hess_upper_accum;

  AD_Func_Opt(
      ad_fun &ag, const std::vector<double> &params_init, bool innerIn = true,
      int num_threads_in = 1,
      const std::vector<std::vector<size_t>> *thread_groups_in = nullptr)
      : AD_Func_Opt(ag, params_init, make_layout(ag, innerIn), innerIn,
                    num_threads_in, thread_groups_in) {}

  int get_nvars() const { return static_cast<int>(nvars_opt); }
  int get_nnz() const { return static_cast<int>(Htemplate.nonZeros()); }
  bool openmp_enabled() const { return num_threads > 1; }
  int openmp_threads() const { return num_threads; }

  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX> &x, double &f) {
    update_gamma_from_x(x);

    double f_sum = 0.0;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      std::vector<double> params_local(parameters.begin(), parameters.end());
      const std::vector<size_t> &shard_group = shard_group_for_thread();

      for (size_t s : shard_group) {
        adlaplace_adpack_handle *h = shards->fun[s];
        GroupPack &pack = *pack_ctx(h->ctx);
        if (!adlaplace_shard_thread_ok(pack)) {
          adlaplace_debug_record_mismatch(
              s, pack.owner_thread,
              static_cast<std::size_t>(omp_get_thread_num()),
              inner ? "inner_opt get_f" : "inner_opt get_f");
          continue;
        }
        (void)h->api->f(h->ctx, params_local.data(), &f_local);
      }

#pragma omp critical(f_sum)
      {
        f_sum -= f_local;
      }
    }

    adlaplace_debug_raise_if_any(inner ? "inner_opt get_f" : "inner_opt get_f");
    f = f_sum;
  }

  template <class DerivedX, class DerivedG>
  void get_fdf(const Eigen::MatrixBase<DerivedX> &x, double &f,
               Eigen::MatrixBase<DerivedG> &g) {
    update_gamma_from_x(x);

    f = 0.0;
    auto &gout = g.derived();
    gout.setZero();

    std::vector<double> grad_full(Nparams, 0.0);
    const bool inner_flag = inner;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      std::vector<double> grad_local(Nparams, 0.0);
      std::vector<double> params_local(parameters.begin(), parameters.end());
      const std::vector<size_t> &shard_group = shard_group_for_thread();

      for (size_t s : shard_group) {
        adlaplace_adpack_handle *h = shards->fun[s];
        GroupPack &pack = *pack_ctx(h->ctx);
        if (!adlaplace_shard_thread_ok(pack)) {
          adlaplace_debug_record_mismatch(
              s, pack.owner_thread,
              static_cast<std::size_t>(omp_get_thread_num()),
              inner ? "inner_opt get_fdf" : "inner_opt get_fdf");
          adlaplace_debug_note_grad_mismatch(grad_local.data(), Nparams);
          continue;
        }
        (void)h->api->f_grad(h->ctx, params_local.data(), &inner_flag, &f_local,
                             grad_local.data());
      }

#pragma omp critical(grad_sum)
      {
        f -= f_local;
        for (size_t k = 0; k < Nparams; ++k) {
          grad_full[k] -= grad_local[k];
        }
      }
    }

    adlaplace_debug_raise_if_any(inner ? "inner_opt get_fdf"
                                       : "inner_opt get_fdf");
    const size_t gsize = static_cast<size_t>(gout.size());
    const size_t ncopy = gsize < nvars_opt ? gsize : nvars_opt;
    for (size_t k = 0; k < ncopy; ++k) {
      gout[static_cast<Eigen::Index>(k)] = grad_full[var_offset + k];
    }
  }

  template <class DerivedX, class DerivedG>
  void get_fdfh(const Eigen::MatrixBase<DerivedX> &x, double &f,
                Eigen::MatrixBase<DerivedG> &g,
                Eigen::SparseMatrix<double> &H) {
    update_gamma_from_x(x);

    f = 0.0;
    auto &gout = g.derived();
    if (gout.size() > 0)
      gout.setZero();

    std::fill(hess_upper_accum.begin(), hess_upper_accum.end(), 0.0);
    std::vector<double> grad_full(Nparams, 0.0);
    const bool inner_flag = inner;
    const char *phase = inner ? "inner_opt get_fdfh" : "outer get_fdfh";
    int api_err_shard = -1;
    int api_err_rc = 0;

#ifdef DEBUG_EXTRA
    std::vector<double> f_local_log(static_cast<size_t>(num_threads), -99.0);
#endif

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      std::vector<double> grad_local(Nparams, 0.0);
      std::vector<double> hess_local(hess_upper_accum.size(), 0.0);
      std::vector<double> params_local(parameters.begin(), parameters.end());
      const std::vector<size_t> &shard_group = shard_group_for_thread();

      for (size_t s : shard_group) {
        adlaplace_adpack_handle *h = shards->fun[s];
        GroupPack &pack = *pack_ctx(h->ctx);
        if (!adlaplace_shard_thread_ok(pack)) {
          adlaplace_debug_record_mismatch(
              s, pack.owner_thread,
              static_cast<std::size_t>(omp_get_thread_num()), phase);
          adlaplace_debug_note_grad_mismatch(grad_local.data(), Nparams);
          continue;
        }
        if (!h->api->f_grad_hess) {
#pragma omp critical(adlaplace_api_err)
          {
            if (api_err_shard < 0) {
              api_err_shard = static_cast<int>(s);
              api_err_rc = -1;
            }
          }
          continue;
        }
        int *map_ptr = const_cast<int *>(hess_maps[s].data());
        std::fill(hess_local.begin(), hess_local.end(), 0.0);
        const int rc = h->api->f_grad_hess(
            h->ctx, params_local.data(), &inner_flag, &f_local,
            grad_local.data(), hess_local.data(), map_ptr);
        if (rc != 0) {
#pragma omp critical(adlaplace_api_err)
          {
            if (api_err_shard < 0) {
              api_err_shard = static_cast<int>(s);
              api_err_rc = rc;
            }
          }
          continue;
        }

#pragma omp critical(hess_sum)
        {
          for (size_t k = 0; k < hess_upper_accum.size(); ++k) {
            hess_upper_accum[k] -= hess_local[k];
          }
        }
      }

#pragma omp critical(grad_sum)
      {
        f -= f_local;
        for (size_t k = 0; k < Nparams; ++k) {
          grad_full[k] -= grad_local[k];
        }
      }

#ifdef DEBUG_EXTRA
      f_local_log[static_cast<size_t>(omp_get_thread_num())] = f_local;
#endif
    }

#ifdef DEBUG_EXTRA
    Rcpp::Rcout << phase << " f_local:\n";
    for (int t = 0; t < num_threads; ++t) {
      Rcpp::Rcout << "  thread=" << t
                  << " f_local=" << f_local_log[static_cast<size_t>(t)] << "\n";
    }
#endif

    if (api_err_shard >= 0) {
      if (api_err_rc < 0) {
        Rcpp::stop("%s: api->f_grad_hess is NULL for shard %d", phase,
                   api_err_shard);
      }
      Rcpp::stop("%s: api->f_grad_hess failed for shard %d with code %d", phase,
                 api_err_shard, api_err_rc);
    }

    adlaplace_debug_raise_if_any(phase);
    if (gout.size() > 0) {
      const size_t gsize = static_cast<size_t>(gout.size());
      const size_t ncopy = gsize < nvars_opt ? gsize : nvars_opt;
      for (size_t k = 0; k < ncopy; ++k) {
        gout[static_cast<Eigen::Index>(k)] = grad_full[var_offset + k];
      }
    }

    if (H.rows() != Htemplate.rows() || H.cols() != Htemplate.cols() ||
        H.nonZeros() != Htemplate.nonZeros()) {
      H = Htemplate.cast<double>();
      H.makeCompressed();
    }

    double *Hx = H.valuePtr();
    const int *cell_id = Htemplate.valuePtr();
    const Eigen::Index nz = H.nonZeros();
    for (Eigen::Index t = 0; t < nz; ++t) {
      Hx[t] = hess_upper_accum[static_cast<size_t>(cell_id[t])];
    }
  }

  template <class DerivedX>
  void get_hessian(const Eigen::MatrixBase<DerivedX> &x,
                   Eigen::SparseMatrix<double> &H) {
    double f_dummy = 0.0;
    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(nvars_opt));
    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedH>
  void get_hessian(const Eigen::MatrixBase<DerivedX> &x,
                   Eigen::SparseMatrixBase<DerivedH> &Hbase) {
    DerivedH &H = Hbase.derived();
    if (H.rows() != Htemplate.rows() || H.cols() != Htemplate.cols() ||
        H.nonZeros() != Htemplate.nonZeros()) {
      H = Htemplate.cast<double>();
      H.makeCompressed();
    }
    double f_dummy = 0.0;
    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(nvars_opt));
    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedG>
  void get_df(const Eigen::MatrixBase<DerivedX> &x,
              Eigen::MatrixBase<DerivedG> &g) {
    double f_dummy = 0.0;
    get_fdf(x, f_dummy, g);
  }

private:
  struct Layout {
    size_t Nparams;
    size_t n_shards_backend;
    size_t nvars_opt;
    size_t var_offset;
    hessian_template Htemplate;
    std::vector<std::vector<int>> hess_maps;
  };

  static std::vector<std::vector<int>> build_hess_maps(ad_fun &ag, bool innerIn,
                                                       size_t n_shards) {

    const hessian_map_view &hv = innerIn ? ag.map_inner : ag.map_outer;
    const R_xlen_t p_len = hv.p.size();
    if (p_len != static_cast<R_xlen_t>(n_shards + 1)) {
      Rcpp::stop("hessian map p length %d but expected n_shards+1=%d",
                 static_cast<int>(p_len), static_cast<int>(n_shards + 1));
    }

    std::vector<std::vector<int>> maps(n_shards);
    for (size_t s = 0; s < n_shards; ++s) {
      adlaplace_adpack_handle *h = ag.fun[s];
      if (!h || !h->api || !h->api->get_sizes) {
        Rcpp::stop("ad_fun.fun[%d] missing get_sizes", static_cast<int>(s));
      }

      int n_inner = 0;
      int n_outer = 0;
      int n_beta = 0;
      int n_theta = 0;
      int nnz_grad_inner = 0;
      int nnz_grad_outer = 0;
      int nnz_hes_inner = 0;
      int nnz_hes_outer = 0;
      if (h->api->get_sizes(h->ctx, &n_inner, &n_outer, &n_beta, &n_theta,
                            &nnz_grad_inner, &nnz_grad_outer, &nnz_hes_inner,
                            &nnz_hes_outer) != 0) {
        Rcpp::stop("get_sizes failed for shard %d", static_cast<int>(s));
      }

      const int nnz_hes = innerIn ? nnz_hes_inner : nnz_hes_outer;
      maps[s].assign(static_cast<size_t>(nnz_hes), -1);

      const int col_start = hv.p[static_cast<R_xlen_t>(s)];
      const int col_end = hv.p[static_cast<R_xlen_t>(s + 1)];
      for (int k = col_start; k < col_end; ++k) {
        const int loc = hv.local[k];
        const int glob = hv.global[k];
        if (loc < 0 || loc >= nnz_hes) {
          Rcpp::stop(
              "hessian map local index %d out of range [0, %d) for shard %d",
              loc, nnz_hes, static_cast<int>(s));
        }
        maps[s][static_cast<size_t>(loc)] = glob;
      }

      for (int loc = 0; loc < nnz_hes; ++loc) {
        if (maps[s][static_cast<size_t>(loc)] < 0) {
          Rcpp::stop("hessian map missing entry for shard %d local index %d",
                     static_cast<int>(s), loc);
        }
      }
    }
    return maps;
  }

  static Layout make_layout(ad_fun &ag, bool innerIn) {
    if (ag.fun.empty())
      Rcpp::stop("ad_fun.fun is empty");

    const size_t Nbeta = static_cast<size_t>(ag.sizes.named("beta"));
    const size_t Ngamma = static_cast<size_t>(ag.sizes.named("gamma"));
    const size_t n_shards = ag.fun.size();
    adlaplace_adpack_handle *h0 = ag.fun[0];
    const size_t Nparams = pack_ctx(h0->ctx)->x.size();

    const size_t nvars_opt = innerIn ? Ngamma : Nparams;
    const size_t var_offset = innerIn ? Nbeta : 0;

    const hessian_template &tpl = innerIn ? ag.hessian_inner : ag.hessian_outer;
    if (tpl.nonZeros() == 0) {
      Rcpp::stop("%s Hessian template missing on ad_fun; attach hessian_map "
                 "result first",
                 innerIn ? "inner" : "outer");
    }
    if (static_cast<size_t>(tpl.rows()) != nvars_opt) {
      Rcpp::stop("%s Hessian nrow=%d but expected nvars=%d",
                 innerIn ? "inner" : "outer", static_cast<int>(tpl.rows()),
                 static_cast<int>(nvars_opt));
    }

    Layout out;
    out.Nparams = Nparams;
    out.n_shards_backend = n_shards;
    out.nvars_opt = nvars_opt;
    out.var_offset = var_offset;
    out.Htemplate = tpl;
    out.hess_maps = build_hess_maps(ag, innerIn, n_shards);
    return out;
  }

  static std::vector<size_t> make_all_shards(size_t n_shards) {
    std::vector<size_t> out(n_shards);
    for (size_t s = 0; s < n_shards; ++s) {
      out[s] = s;
    }
    return out;
  }

  AD_Func_Opt(ad_fun &ag, const std::vector<double> &params_init, Layout layout,
              bool innerIn, int num_threads_in,
              const std::vector<std::vector<size_t>> *thread_groups_in)
      : shards(&ag), thread_groups(thread_groups_in), inner(innerIn),
        num_threads(num_threads_in > 0 ? num_threads_in : 1),
        Nparams(layout.Nparams), n_shards_backend(layout.n_shards_backend),
        Nbeta_backend(static_cast<size_t>(ag.sizes.named("beta"))),
        Ngamma_backend(static_cast<size_t>(ag.sizes.named("gamma"))),
        Ntheta_backend(static_cast<size_t>(ag.sizes.named("theta"))),
        nvars_opt(layout.nvars_opt), var_offset(layout.var_offset),
        Htemplate(layout.Htemplate), hess_maps(std::move(layout.hess_maps)),
        all_shards(make_all_shards(layout.n_shards_backend)), parameters(),
        hess_upper_accum(static_cast<size_t>(Htemplate.nonZeros()), 0.0) {
    if (thread_groups &&
        static_cast<int>(thread_groups->size()) != num_threads) {
      Rcpp::stop("thread affinity group count %d does not match num_threads %d",
                 static_cast<int>(thread_groups->size()), num_threads);
    }
    parameters.resize(Nparams, 0.0);
    for (size_t j = 0; j < params_init.size() && j < Nparams; ++j) {
      parameters[j] = params_init[j];
    }
  }

  inline const std::vector<size_t> &shard_group_for_thread() const {
    if (!thread_groups) {
      return all_shards;
    }
    const int tid = omp_get_thread_num();
    if (tid < 0 || tid >= static_cast<int>(thread_groups->size())) {
      Rcpp::stop(
          "OpenMP thread %d out of range for thread affinity groups (%d)", tid,
          static_cast<int>(thread_groups->size()));
    }
    return (*thread_groups)[static_cast<size_t>(tid)];
  }

  template <class DerivedX>
  inline void update_gamma_from_x(const Eigen::MatrixBase<DerivedX> &x) {
    if (static_cast<size_t>(x.size()) != nvars_opt) {
      return;
    }
    if (inner) {
      for (size_t k = 0; k < Ngamma_backend; ++k) {
        parameters[Nbeta_backend + k] = x[static_cast<Eigen::Index>(k)];
      }
    } else {
      for (size_t k = 0; k < Nparams; ++k) {
        parameters[k] = x[static_cast<Eigen::Index>(k)];
      }
    }
  }
};

#endif
