#ifndef ADLAPLACE_AD_FUNC_OPT_HPP
#define ADLAPLACE_AD_FUNC_OPT_HPP

#include <Rcpp.h>
#include <Rinternals.h>
#include <cstdlib>
#include <cstdio>

#include <Eigen/Sparse>
#include <numeric>
#include <vector>

#include <omp.h>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/ompad.hpp"
#include "adlaplace/runtime/backend.hpp"
#include "adlaplace/runtime/rviews.hpp"

struct AD_Func_Opt {
  SEXP adFun;
  adlaplace_adpack_handle* handle;
  bool inner;
  int num_threads;
  bool cppad_setup_done;
  std::vector<std::vector<GroupPack>> adfun_local;
  std::vector<BackendContext> ctx_local;

  size_t Nparams;
  size_t Ngroups_backend;
  size_t Nbeta_backend;
  size_t Ngamma_backend;
  size_t Ntheta_backend;
  size_t nvars_opt;
  size_t var_offset;

  // maps full symmetric Hessian nz -> index in backend upper-triangle buffer
  Eigen::SparseMatrix<int, Eigen::ColMajor, int> Htemplate;
  std::vector<int> h_index_upper;

  // full parameter vector reused across calls
  std::vector<double> parameters;
  std::vector<double> hess_upper_accum;

  AD_Func_Opt(
    SEXP adPackIn,
    const std::vector<double>& params_init,
    bool innerIn = true,
    int num_threads_in = 1
    )
  : adFun(adPackIn),
  handle(get_handle(adPackIn)),
  inner(innerIn),
  num_threads(num_threads_in > 0 ? num_threads_in : 1),
  cppad_setup_done(false),
  Nparams(0),
  Ngroups_backend(0),
  Nbeta_backend(0),
  Ngamma_backend(0),
  Ntheta_backend(0),
  nvars_opt(0),
  var_offset(0)
  {
    BackendContext* base_ctx = static_cast<BackendContext*>(handle->ctx);
    if (!base_ctx || !base_ctx->adFun) {
      Rcpp::stop("adFun backend context is missing");
    }
    rebuild_thread_contexts_from_base(*base_ctx);
    void* ctx0 = static_cast<void*>(&ctx_local[0]);

    const int rc_sizes = handle->api->get_sizes(
      ctx0, &Nparams, &Ngroups_backend, &Nbeta_backend, &Ngamma_backend, &Ntheta_backend
      );
    if (rc_sizes != 0) {
      Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
    }
    parameters.resize(Nparams, 0.0);
    for (size_t j = 0; j < params_init.size() && j < Nparams; ++j) {
      parameters[j] = params_init[j];
    }

    nvars_opt = inner ? Ngamma_backend : Nparams;
    var_offset = inner ? Nbeta_backend : 0;

    const bool inner_flag = inner;
    const int* p = NULL;
    const int* i = NULL;
    size_t p_len = 0;
    size_t i_len = 0;
    const int rc_hess = handle->api->get_hessian(ctx0, &inner_flag, &p, &p_len, &i, &i_len);
    if (rc_hess != 0) {
      Rcpp::stop("backend api->get_hessian failed with code %d", rc_hess);
    }
    if (p_len != (nvars_opt + 1)) {
      Rcpp::stop(
        "%s Hessian p_len=%d but expected nvars+1=%d",
        inner ? "inner" : "outer",
        (int)p_len,
        (int)(nvars_opt + 1)
        );
    }

    hess_upper_accum.resize(i_len, 0.0);
    h_index_upper.resize(i_len);
    std::iota(h_index_upper.begin(), h_index_upper.end(), 0);

    using SpMatI = Eigen::SparseMatrix<int, Eigen::ColMajor, int>;
    Eigen::Map<const SpMatI> Hupper(
      static_cast<int>(nvars_opt),
      static_cast<int>(nvars_opt),
      static_cast<int>(i_len),
      p, i, h_index_upper.data()
      );
    Htemplate = Hupper.selfadjointView<Eigen::Upper>();
    Htemplate.makeCompressed();
  }

  int get_nvars() const { return static_cast<int>(nvars_opt); }
  int get_nnz() const { return static_cast<int>(Htemplate.nonZeros()); }
  bool openmp_enabled() const { return num_threads > 1; }
  int openmp_threads() const { return num_threads; }

  void reset_thread_contexts() {
    BackendContext* base_ctx = static_cast<BackendContext*>(handle->ctx);
    if (!base_ctx || !base_ctx->adFun) {
      return;
    }
    rebuild_thread_contexts_from_base(*base_ctx);
  }



  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX>& x, double& f) {
    reset_thread_contexts();
    ensure_cppad_parallel_setup();

    update_gamma_from_x(x);

    double f_sum = 0.0;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      int gend = static_cast<int>(Ngroups_backend);
      std::vector<double> params_local(parameters.begin(), parameters.end());
      const int tid = omp_get_thread_num();
      const size_t tidx = (static_cast<size_t>(tid) < ctx_local.size()) ? static_cast<size_t>(tid) : 0;
      void* ctx_ptr = static_cast<void*>(&ctx_local[tidx]);

#pragma omp for schedule(static,1)
      for (int g = 0; g < gend; ++g) {
        int gi = g;
        (void)handle->api->f(ctx_ptr, &gi, params_local.data(), &f_local);
        } // g loop
      #pragma omp critical
      {
        f_sum += f_local;
      } // omp critical
    } // pragma omp parallel
    // trustOptim minimizes: return minus log-likelihood
    f = -f_sum;
  }

  template <class DerivedX, class DerivedG>
    void get_fdf(const Eigen::MatrixBase<DerivedX>& x, double& f, Eigen::MatrixBase<DerivedG>& g) {
      reset_thread_contexts();
      ensure_cppad_parallel_setup();

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
        int gend = static_cast<int>(Ngroups_backend);
        std::vector<double> params_local(parameters.begin(), parameters.end());
        const int tid = omp_get_thread_num();
        const size_t tidx = (static_cast<size_t>(tid) < ctx_local.size()) ? static_cast<size_t>(tid) : 0;
        void* ctx_ptr = static_cast<void*>(&ctx_local[tidx]);

#pragma omp for schedule(static,1)
        for (int g = 0; g < gend; ++g) {
          int gi = g;
          (void)handle->api->f_grad(
            ctx_ptr, &gi, params_local.data(), &inner_flag, &f_local, grad_local.data()
            );
        }

        #pragma omp critical
        {
          f += f_local;
          for (size_t k = 0; k < Nparams; ++k) {
            grad_full[k] += grad_local[k];
          }
        }
      }

      f = -f;
      const size_t gsize = static_cast<size_t>(gout.size());
      const size_t ncopy = gsize < nvars_opt ? gsize : nvars_opt;
      for (size_t k = 0; k < ncopy; ++k) {
        gout[static_cast<Eigen::Index>(k)] = -grad_full[var_offset + k];
      }
    }

  template <class DerivedX, class DerivedG>
    void get_fdfh(
      const Eigen::MatrixBase<DerivedX>& x,
      double& f,
      Eigen::MatrixBase<DerivedG>& g,
      Eigen::SparseMatrix<double>& H
      ) {
      reset_thread_contexts();
      ensure_cppad_parallel_setup();

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
        int gend = static_cast<int>(Ngroups_backend);
        std::vector<double> params_local(parameters.begin(), parameters.end());
        const int tid = omp_get_thread_num();
        const size_t tidx = (static_cast<size_t>(tid) < ctx_local.size()) ? static_cast<size_t>(tid) : 0;
        void* ctx_ptr = static_cast<void*>(&ctx_local[tidx]);

        #pragma omp for schedule(static,1)
        for (int g = 0; g < gend; ++g) {
          int gi = g;
          (void)handle->api->f_grad_hess(
            ctx_ptr, &gi, params_local.data(), &inner_flag,
            &f_local, grad_local.data(), hess_local.data()
            );
        }

        #pragma omp critical
        {
          f += f_local;
          for (size_t k = 0; k < Nparams; ++k) {
            grad_full[k] += grad_local[k];
          }
          for (size_t k = 0; k < hess_upper_accum.size(); ++k) {
            hess_upper_accum[k] += hess_local[k];
          }
        } // critical
      } // parallel

      f = -f;
      if (gout.size() > 0) {
        const size_t gsize = static_cast<size_t>(gout.size());
        const size_t ncopy = gsize < nvars_opt ? gsize : nvars_opt;
        for (size_t k = 0; k < ncopy; ++k) {
          gout[static_cast<Eigen::Index>(k)] = -grad_full[var_offset + k];
        }
      }

      if (H.rows() != Htemplate.rows() || H.cols() != Htemplate.cols() || H.nonZeros() != Htemplate.nonZeros()) {
        H = Htemplate.cast<double>();
        H.makeCompressed();
      }

      double* Hx = H.valuePtr();
      const int* map = Htemplate.valuePtr();
      const Eigen::Index nz = H.nonZeros();
      for (Eigen::Index t = 0; t < nz; ++t) {
        Hx[t] = -hess_upper_accum[static_cast<size_t>(map[t])];
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
    inline void ensure_cppad_parallel_setup() {
      if (!cppad_setup_done) {
        cppad_parallel_setup(static_cast<std::size_t>(num_threads));
        cppad_setup_done = true;
      }
    }

    void rebuild_thread_contexts_from_base(const BackendContext& base_ctx) {
      adfun_local.resize(static_cast<size_t>(num_threads));
      ctx_local.resize(static_cast<size_t>(num_threads));
      for (int t = 0; t < num_threads; ++t) {
        adfun_local[static_cast<size_t>(t)] = clone_adfun(*(base_ctx.adFun));
        BackendContext& ctx_t = ctx_local[static_cast<size_t>(t)];
        ctx_t.adFun = &adfun_local[static_cast<size_t>(t)];
        ctx_t.Nparams = base_ctx.Nparams;
        ctx_t.Ngroups = base_ctx.Ngroups;
        ctx_t.Nbeta = base_ctx.Nbeta;
        ctx_t.Ngamma = base_ctx.Ngamma;
        ctx_t.Ntheta = base_ctx.Ntheta;
        ctx_t.hessian_inner = base_ctx.hessian_inner;
        ctx_t.hessian_outer = base_ctx.hessian_outer;
      }
    }

    static GroupPack clone_group_pack(const GroupPack& src) {
      GroupPack dst;

      // Deep-copy AD tape/state and persistent sparsity structures.
      dst.fun = src.fun;
      dst.pattern_grad = src.pattern_grad;
      dst.pattern_grad_inner = src.pattern_grad_inner;
      dst.pattern_hessian = src.pattern_hessian;
      dst.pattern_hessian_inner = src.pattern_hessian_inner;
      dst.w = src.w;
      dst.wthree = src.wthree;
      dst.direction_zeros = src.direction_zeros;
      dst.direction = src.direction;
      dst.unused_pattern = src.unused_pattern;
      dst.x = src.x;

      // Keep work caches thread-local and freshly initialized.
      dst.work_grad = CppAD::sparse_jac_work();
      dst.work_hess = CppAD::sparse_hes_work();
      dst.work_inner_grad = CppAD::sparse_jac_work();
      dst.work_inner_hess = CppAD::sparse_hes_work();

      return dst;
    }

    static std::vector<GroupPack> clone_adfun(const std::vector<GroupPack>& src) {
      std::vector<GroupPack> out;
      out.reserve(src.size());
      for (const GroupPack& gp : src) {
        out.push_back(clone_group_pack(gp));
      }
      return out;
    }

    static adlaplace_adpack_handle* get_handle(SEXP handle_sexp) {
      SEXP handle_ptr = handle_sexp;
      if (TYPEOF(handle_sexp) == VECSXP) {
        Rcpp::List maybe_list(handle_sexp);
        if (maybe_list.containsElementNamed("adFun")) {
          handle_ptr = maybe_list["adFun"];
        }
      }

      adlaplace_adpack_handle* h =
      static_cast<adlaplace_adpack_handle*>(R_ExternalPtrAddr(handle_ptr));
      if (!h) Rcpp::stop("adFun handle is NULL");
      if (!h->api) Rcpp::stop("adFun handle api is NULL");
      if (!h->ctx) Rcpp::stop("adFun handle ctx is NULL");
      if (!h->api->f || !h->api->f_grad || !h->api->f_grad_hess ||
        !h->api->get_hessian || !h->api->get_sizes || !h->api->trace_hinv_t) {
        Rcpp::stop("adFun handle has incomplete API function table");
    }
    return h;
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
