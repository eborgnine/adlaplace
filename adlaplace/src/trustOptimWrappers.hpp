#ifndef ADLAPLACE_AD_FUNC_OPT_HPP
#define ADLAPLACE_AD_FUNC_OPT_HPP

#include <Rcpp.h>
#include <Rinternals.h>

#include <Eigen/Sparse>
#include <numeric>
#include <vector>

#include <omp.h>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/runtime/interfaces.hpp"
#include "adlaplace/runtime/rviews.hpp"

struct AD_Func_Opt {
  ad_groups* const groups;
  const bool inner;
  const int num_threads;

  const size_t Nparams;
  const size_t Ngroups_backend;
  const size_t Nbeta_backend;
  const size_t Ngamma_backend;
  const size_t Ntheta_backend;
  const size_t nvars_opt;
  const size_t var_offset;

  const Eigen::SparseMatrix<int, Eigen::ColMajor, int> Htemplate;
  const std::vector<std::vector<int>> hess_maps;

  std::vector<double> parameters;
  std::vector<double> hess_upper_accum;

  AD_Func_Opt(
    ad_groups& ag,
    const std::vector<double>& params_init,
    bool innerIn = true,
    int num_threads_in = 1
    )
  : AD_Func_Opt(ag, params_init, make_layout(ag, innerIn), innerIn, num_threads_in)
  {}

  int get_nvars() const { return static_cast<int>(nvars_opt); }
  int get_nnz() const { return static_cast<int>(Htemplate.nonZeros()); }
  bool openmp_enabled() const { return num_threads > 1; }
  int openmp_threads() const { return num_threads; }

  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX>& x, double& f) {
    update_gamma_from_x(x);

    double f_sum = 0.0;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      const int gend = static_cast<int>(Ngroups_backend);
      std::vector<double> params_local(parameters.begin(), parameters.end());

#pragma omp for schedule(static,1) nowait
      for (int g = 0; g < gend; ++g) {
        adlaplace_adpack_handle* h = groups->fun[static_cast<size_t>(g)];
        (void)h->api->f(h->ctx, params_local.data(), &f_local);
      }

#pragma omp critical
      {
        f_sum -= f_local;
      }
    }

    f = f_sum;
  }

  template <class DerivedX, class DerivedG>
  void get_fdf(const Eigen::MatrixBase<DerivedX>& x, double& f, Eigen::MatrixBase<DerivedG>& g) {
    update_gamma_from_x(x);

    f = 0.0;
    auto& gout = g.derived();
    gout.setZero();

    std::vector<double> grad_full(Nparams, 0.0);
    const bool inner_flag = inner;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      std::vector<double> grad_local(Nparams, 0.0);
      const int gend = static_cast<int>(Ngroups_backend);
      std::vector<double> params_local(parameters.begin(), parameters.end());

#pragma omp for schedule(static,1) nowait
      for (int g = 0; g < gend; ++g) {
        adlaplace_adpack_handle* h = groups->fun[static_cast<size_t>(g)];
        (void)h->api->f_grad(
          h->ctx, params_local.data(), &inner_flag, &f_local, grad_local.data()
        );
      }

#pragma omp critical
      {
        f -= f_local;
        for (size_t k = 0; k < Nparams; ++k) {
          grad_full[k] -= grad_local[k];
        }
      }
    }

    const size_t gsize = static_cast<size_t>(gout.size());
    const size_t ncopy = gsize < nvars_opt ? gsize : nvars_opt;
    for (size_t k = 0; k < ncopy; ++k) {
      gout[static_cast<Eigen::Index>(k)] = grad_full[var_offset + k];
    }
  }

  template <class DerivedX, class DerivedG>
  void get_fdfh(
    const Eigen::MatrixBase<DerivedX>& x,
    double& f,
    Eigen::MatrixBase<DerivedG>& g,
    Eigen::SparseMatrix<double>& H
  ) {
    update_gamma_from_x(x);

    f = 0.0;
    auto& gout = g.derived();
    if (gout.size() > 0) gout.setZero();

    std::fill(hess_upper_accum.begin(), hess_upper_accum.end(), 0.0);
    std::vector<double> grad_full(Nparams, 0.0);
    const bool inner_flag = inner;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      std::vector<double> grad_local(Nparams, 0.0);
      std::vector<double> hess_local(hess_upper_accum.size(), 0.0);
      const int gend = static_cast<int>(Ngroups_backend);
      std::vector<double> params_local(parameters.begin(), parameters.end());

#if ADLAPLACE_HAVE_KMPC_DISPATCH_DEINIT
#pragma omp for schedule(dynamic) nowait
#else
#pragma omp for schedule(static,1) nowait
#endif
      for (int g = 0; g < gend; ++g) {
        adlaplace_adpack_handle* h = groups->fun[static_cast<size_t>(g)];
        int* map_ptr = const_cast<int*>(hess_maps[static_cast<size_t>(g)].data());
        std::fill(hess_local.begin(), hess_local.end(), 0.0);
        (void)h->api->f_grad_hess(
          h->ctx, params_local.data(), &inner_flag,
          &f_local, grad_local.data(), hess_local.data(), map_ptr
        );

#pragma omp critical
        {
          for (size_t k = 0; k < hess_upper_accum.size(); ++k) {
            hess_upper_accum[k] -= hess_local[k];
          }
        }
      }

#pragma omp critical
      {
        f -= f_local;
        for (size_t k = 0; k < Nparams; ++k) {
          grad_full[k] -= grad_local[k];
        }
      }
    }

    if (gout.size() > 0) {
      const size_t gsize = static_cast<size_t>(gout.size());
      const size_t ncopy = gsize < nvars_opt ? gsize : nvars_opt;
      for (size_t k = 0; k < ncopy; ++k) {
        gout[static_cast<Eigen::Index>(k)] = grad_full[var_offset + k];
      }
    }

    if (H.rows() != Htemplate.rows() || H.cols() != Htemplate.cols() || H.nonZeros() != Htemplate.nonZeros()) {
      H = Htemplate.cast<double>();
      H.makeCompressed();
    }

    double* Hx = H.valuePtr();
    const int* cell_id = Htemplate.valuePtr();
    const Eigen::Index nz = H.nonZeros();
    for (Eigen::Index t = 0; t < nz; ++t) {
      Hx[t] = hess_upper_accum[static_cast<size_t>(cell_id[t])];
    }
  }

  template <class DerivedX>
  void get_hessian(const Eigen::MatrixBase<DerivedX>& x, Eigen::SparseMatrix<double>& H) {
    double f_dummy = 0.0;
    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(nvars_opt));
    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedH>
  void get_hessian(const Eigen::MatrixBase<DerivedX>& x, Eigen::SparseMatrixBase<DerivedH>& Hbase) {
    DerivedH& H = Hbase.derived();
    if (H.rows() != Htemplate.rows() || H.cols() != Htemplate.cols() || H.nonZeros() != Htemplate.nonZeros()) {
      H = Htemplate.cast<double>();
      H.makeCompressed();
    }
    double f_dummy = 0.0;
    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(nvars_opt));
    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedG>
  void get_df(const Eigen::MatrixBase<DerivedX>& x, Eigen::MatrixBase<DerivedG>& g) {
    double f_dummy = 0.0;
    get_fdf(x, f_dummy, g);
  }

private:
  struct Layout {
    size_t Nparams;
    size_t Ngroups_backend;
    size_t nvars_opt;
    size_t var_offset;
    hessian_template Htemplate;
    std::vector<std::vector<int>> hess_maps;
  };

  static std::vector<std::vector<int>> build_hess_maps(
    ad_groups& ag,
    bool innerIn,
    size_t Ngroups) {

    const hessian_map_view& hv = innerIn ? ag.map_inner : ag.map_outer;
    const R_xlen_t p_len = hv.p.size();
    if (p_len != static_cast<R_xlen_t>(Ngroups + 1)) {
      Rcpp::stop(
        "hessian map p length %d but expected Ngroups+1=%d",
        static_cast<int>(p_len),
        static_cast<int>(Ngroups + 1)
      );
    }

    std::vector<std::vector<int>> maps(Ngroups);
    for (size_t g = 0; g < Ngroups; ++g) {
      adlaplace_adpack_handle* h = ag.fun[g];
      if (!h || !h->api || !h->api->get_sparse_sizes) {
        Rcpp::stop("ad_groups.fun[%d] missing get_sparse_sizes", static_cast<int>(g));
      }

      int n_inner = 0;
      int n_outer = 0;
      int nnz_grad_inner = 0;
      int nnz_grad_outer = 0;
      int nnz_hes_inner = 0;
      int nnz_hes_outer = 0;
      if (h->api->get_sparse_sizes(
            h->ctx, &n_inner, &n_outer,
            &nnz_grad_inner, &nnz_grad_outer,
            &nnz_hes_inner, &nnz_hes_outer) != 0) {
        Rcpp::stop("get_sparse_sizes failed for group %d", static_cast<int>(g));
      }

      const int nnz_hes = innerIn ? nnz_hes_inner : nnz_hes_outer;
      maps[g].assign(static_cast<size_t>(nnz_hes), -1);

      const int col_start = hv.p[static_cast<R_xlen_t>(g)];
      const int col_end = hv.p[static_cast<R_xlen_t>(g + 1)];
      for (int k = col_start; k < col_end; ++k) {
        const int loc = hv.local[k];
        const int glob = hv.global[k];
        if (loc < 0 || loc >= nnz_hes) {
          Rcpp::stop(
            "hessian map local index %d out of range [0, %d) for group %d",
            loc, nnz_hes, static_cast<int>(g)
          );
        }
        maps[g][static_cast<size_t>(loc)] = glob;
      }

      for (int loc = 0; loc < nnz_hes; ++loc) {
        if (maps[g][static_cast<size_t>(loc)] < 0) {
          Rcpp::stop(
            "hessian map missing entry for group %d local index %d",
            static_cast<int>(g),
            loc
          );
        }
      }
    }
    return maps;
  }

  static Layout make_layout(ad_groups& ag, bool innerIn) {
    if (ag.fun.empty()) Rcpp::stop("ad_groups.fun is empty");

    const size_t Nbeta = static_cast<size_t>(ag.sizes.named("beta"));
    const size_t Ngamma = static_cast<size_t>(ag.sizes.named("gamma"));
    const size_t Ngroups = ag.fun.size();
    adlaplace_adpack_handle* h0 = ag.fun[0];
    const size_t Nparams = pack_ctx(h0->ctx)->x.size();

    const size_t nvars_opt = innerIn ? Ngamma : Nparams;
    const size_t var_offset = innerIn ? Nbeta : 0;

    const hessian_template& tpl = innerIn ? ag.hessian_inner : ag.hessian_outer;
    if (tpl.nonZeros() == 0) {
      Rcpp::stop(
        "%s Hessian template missing on ad_groups; attach hessian_map result first",
        innerIn ? "inner" : "outer"
      );
    }
    if (static_cast<size_t>(tpl.rows()) != nvars_opt) {
      Rcpp::stop(
        "%s Hessian nrow=%d but expected nvars=%d",
        innerIn ? "inner" : "outer",
        static_cast<int>(tpl.rows()),
        static_cast<int>(nvars_opt)
      );
    }

    Layout out;
    out.Nparams = Nparams;
    out.Ngroups_backend = Ngroups;
    out.nvars_opt = nvars_opt;
    out.var_offset = var_offset;
    out.Htemplate = tpl;
    out.hess_maps = build_hess_maps(ag, innerIn, Ngroups);
    return out;
  }

  AD_Func_Opt(
    ad_groups& ag,
    const std::vector<double>& params_init,
    Layout layout,
    bool innerIn,
    int num_threads_in
    )
  : groups(&ag),
    inner(innerIn),
    num_threads(num_threads_in > 0 ? num_threads_in : 1),
    Nparams(layout.Nparams),
    Ngroups_backend(layout.Ngroups_backend),
    Nbeta_backend(static_cast<size_t>(ag.sizes.named("beta"))),
    Ngamma_backend(static_cast<size_t>(ag.sizes.named("gamma"))),
    Ntheta_backend(static_cast<size_t>(ag.sizes.named("theta"))),
    nvars_opt(layout.nvars_opt),
    var_offset(layout.var_offset),
    Htemplate(layout.Htemplate),
    hess_maps(std::move(layout.hess_maps)),
    parameters(),
    hess_upper_accum(static_cast<size_t>(Htemplate.nonZeros()), 0.0)
  {
    parameters.resize(Nparams, 0.0);
    for (size_t j = 0; j < params_init.size() && j < Nparams; ++j) {
      parameters[j] = params_init[j];
    }
  }

  template <class DerivedX>
  inline void update_gamma_from_x(const Eigen::MatrixBase<DerivedX>& x) {
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
