#ifndef ADLAPLACE_AD_FUNC_OPT_HPP
#define ADLAPLACE_AD_FUNC_OPT_HPP

#include <Rcpp.h>
#include <Rinternals.h>

#include <Eigen/Sparse>
#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "adlaplace/adpack_handle.h"
#include "adlaplace/rviews.hpp"

struct AD_Func_Opt {
  SEXP adPack;
  adlaplace_adpack_handle* handle;

  size_t Nparams;
  size_t Ngroups_backend;
  size_t Nbeta_backend;
  size_t Ngamma_backend;
  size_t Ntheta_backend;

  // maps full symmetric Hessian nz -> index in backend upper-triangle buffer
  Eigen::SparseMatrix<int, Eigen::ColMajor, int> Htemplate;
  std::vector<int> h_index_upper;

  // full parameter vector reused across calls
  std::vector<double> parameters;
  std::vector<double> hess_upper_accum;

  AD_Func_Opt(SEXP adPackIn, const std::vector<double>& params_init)
    : adPack(adPackIn),
      handle(get_handle(adPackIn)),
      Nparams(0),
      Ngroups_backend(0),
      Nbeta_backend(0),
      Ngamma_backend(0),
      Ntheta_backend(0)
  {
    const int rc_sizes = handle->api->get_sizes(
      handle->ctx, &Nparams, &Ngroups_backend, &Nbeta_backend, &Ngamma_backend, &Ntheta_backend
    );
    if (rc_sizes != 0) {
      Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
    }
    parameters.resize(Nparams, 0.0);
    for (size_t j = 0; j < params_init.size() && j < Nparams; ++j) {
      parameters[j] = params_init[j];
    }

    const bool inner = true;
    const int* p = NULL;
    const int* i = NULL;
    size_t p_len = 0;
    size_t i_len = 0;
    const int rc_hess = handle->api->get_hessian(handle->ctx, &inner, &p, &p_len, &i, &i_len);
    if (rc_hess != 0) {
      Rcpp::stop("backend api->get_hessian failed with code %d", rc_hess);
    }
    if (p_len != (Ngamma_backend + 1)) {
      Rcpp::stop("inner Hessian p_len=%d but expected Ngamma+1=%d", (int)p_len, (int)(Ngamma_backend + 1));
    }

    hess_upper_accum.resize(i_len, 0.0);
    h_index_upper.resize(i_len);
    std::iota(h_index_upper.begin(), h_index_upper.end(), 0);

    using SpMatI = Eigen::SparseMatrix<int, Eigen::ColMajor, int>;
    Eigen::Map<const SpMatI> Hupper(
      static_cast<int>(Ngamma_backend),
      static_cast<int>(Ngamma_backend),
      static_cast<int>(i_len),
      p, i, h_index_upper.data()
    );
    Htemplate = Hupper.selfadjointView<Eigen::Upper>();
    Htemplate.makeCompressed();
  }

  int get_nvars() const { return static_cast<int>(Ngamma_backend); }
  int get_nnz() const { return static_cast<int>(Htemplate.nonZeros()); }

  template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX>& x, double& f) {
    update_gamma_from_x(x);

    double f_sum = 0.0;
    for (size_t g = 0; g < Ngroups_backend; ++g) {
      int gi = static_cast<int>(g);
      int rc = handle->api->f(handle->ctx, &gi, parameters.data(), &f_sum);
      if (rc != 0) {
        Rcpp::stop("backend api->f failed for group %d with code %d", gi, rc);
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
    const bool inner = true;

    for (size_t g = 0; g < Ngroups_backend; ++g) {
      int gi = static_cast<int>(g);
      int rc = handle->api->f_grad(
        handle->ctx, &gi, parameters.data(), &inner, &f, grad_full.data()
      );
      if (rc != 0) {
        Rcpp::stop("backend api->f_grad failed for group %d with code %d", gi, rc);
      }
    }

    f = -f;
    for (size_t k = 0; k < Ngamma_backend; ++k) {
      gout[static_cast<Eigen::Index>(k)] = -grad_full[Nbeta_backend + k];
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
    const bool inner = true;

    for (size_t g = 0; g < Ngroups_backend; ++g) {
      int gi = static_cast<int>(g);
      int rc = handle->api->f_grad_hess(
        handle->ctx, &gi, parameters.data(), &inner, &f, grad_full.data(), hess_upper_accum.data()
      );
      if (rc != 0) {
        Rcpp::stop("backend api->f_grad_hess failed for group %d with code %d", gi, rc);
      }
    }

    f = -f;
    if (gout.size() > 0) {
      if (static_cast<size_t>(gout.size()) != Ngamma_backend) {
        Rcpp::stop("gradient buffer has size %d but expected %d", (int)gout.size(), (int)Ngamma_backend);
      }
      for (size_t k = 0; k < Ngamma_backend; ++k) {
        gout[static_cast<Eigen::Index>(k)] = -grad_full[Nbeta_backend + k];
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
    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(Ngamma_backend));
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
    Eigen::VectorXd g_dummy(static_cast<Eigen::Index>(Ngamma_backend));
    get_fdfh(x, f_dummy, g_dummy, H);
  }

  template <class DerivedX, class DerivedG>
  void get_df(const Eigen::MatrixBase<DerivedX>& x, Eigen::MatrixBase<DerivedG>& g) {
    double f_dummy = 0.0;
    get_fdf(x, f_dummy, g);
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
    if (!h) Rcpp::stop("adPack handle is NULL");
    if (!h->api) Rcpp::stop("adPack handle api is NULL");
    if (!h->ctx) Rcpp::stop("adPack handle ctx is NULL");
    if (!h->api->f || !h->api->f_grad || !h->api->f_grad_hess || !h->api->get_hessian || !h->api->get_sizes) {
      Rcpp::stop("adPack handle has incomplete API function table");
    }
    return h;
  }

  template <class DerivedX>
  inline void update_gamma_from_x(const Eigen::MatrixBase<DerivedX>& x) {
    for (size_t k = 0; k < Ngamma_backend; ++k) {
      parameters[Nbeta_backend + k] = x[static_cast<Eigen::Index>(k)];
    }
  }
};

#endif
