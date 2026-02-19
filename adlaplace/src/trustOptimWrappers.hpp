#ifndef ADLAPLACE_AD_FUNC_OPT_HPP
#define ADLAPLACE_AD_FUNC_OPT_HPP

#include <Rcpp.h>
#include <Rinternals.h>

#include <Eigen/Sparse>
#include <numeric>
#include <vector>

#include <omp.h>

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/runtime/rviews.hpp"

struct AD_Func_Opt {
  SEXP adFun;
  adlaplace_adpack_handle* handle;
  bool inner;
  int num_threads;

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
    Nparams(0),
    Ngroups_backend(0),
    Nbeta_backend(0),
    Ngamma_backend(0),
    Ntheta_backend(0),
    nvars_opt(0),
    var_offset(0)
  {
    void* ctx = handle->ctx;

    const int rc_sizes = handle->api->get_sizes(
      ctx, &Nparams, &Ngroups_backend, &Nbeta_backend, &Ngamma_backend, &Ntheta_backend
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
    const int rc_hess = handle->api->get_hessian(ctx, &inner_flag, &p, &p_len, &i, &i_len);
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

  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX>& x, double& f) {
    update_gamma_from_x(x);

    double f_sum = 0.0;

#pragma omp parallel num_threads(num_threads)
    {
      double f_local = 0.0;
      const int gend = static_cast<int>(Ngroups_backend);
      std::vector<double> params_local(parameters.begin(), parameters.end());

#pragma omp for schedule(static,1) 
      for (int g = 0; g < gend; ++g) {
        int gi = g;
        (void)handle->api->f(handle->ctx, &gi, params_local.data(), &f_local);
      }

#pragma omp critical
      {
        f_sum += f_local;
      }
    }

    // trustOptim minimizes: return minus log-likelihood
    f = -f_sum;
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

#pragma omp for schedule(static,1) 
      for (int g = 0; g < gend; ++g) {
        int gi = g;
        (void)handle->api->f_grad(
          handle->ctx, &gi, params_local.data(), &inner_flag, &f_local, grad_local.data()
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
#pragma omp for schedule(dynamic) 
#else      
#pragma omp for schedule(static,1) 
#endif      
      for (int g = 0; g < gend; ++g) {
        int gi = g;
        (void)handle->api->f_grad_hess(
          handle->ctx, &gi, params_local.data(), &inner_flag,
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
      }
    }

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


// empty functions, called by get_current_state
  template <class DerivedX>
  void get_hessian(const Eigen::MatrixBase<DerivedX>& x, Eigen::SparseMatrix<double>& H) {
//    double f_dummy = 0.0;
//    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(nvars_opt));
//    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedH>
  void get_hessian(const Eigen::MatrixBase<DerivedX>& x, Eigen::SparseMatrixBase<DerivedH>& Hbase) {
//    DerivedH& H = Hbase.derived();
//    if (H.rows() != Htemplate.rows() || H.cols() != Htemplate.cols() || H.nonZeros() != Htemplate.nonZeros()) {
//      H = Htemplate.cast<double>();
//      H.makeCompressed();
//    }
//    double f_dummy = 0.0;
//    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(nvars_opt));
//    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedG>
  void get_df(const Eigen::MatrixBase<DerivedX>& x, Eigen::MatrixBase<DerivedG>& g) {
//    double f_dummy = 0.0;
//    get_fdf(x, f_dummy, g);
  }

private:
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
