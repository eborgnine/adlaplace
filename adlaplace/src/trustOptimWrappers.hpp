#ifndef ADLAPLACE_AD_FUNC_OPT_HPP
#define ADLAPLACE_AD_FUNC_OPT_HPP

#include <Rcpp.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <Eigen/Sparse>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "adlaplace/utils.hpp"
#include "adlaplace/adlaplace_api.hpp"

struct AD_Func_Opt {
  SEXP adPack;
  Config config;

// “index templates”
  Eigen::SparseMatrix<int, Eigen::ColMajor, int> Htemplate;

// buffers reused across calls
  CPPAD_TESTVECTOR(double) parameters;
  CPPAD_TESTVECTOR(double) result_hessian;

  std::vector<int> h_index_upper;


  AD_Func_Opt(SEXP adPackIn, const Config& configIn)
  : adPack(adPackIn),
  config(configIn),
  parameters(config.Nparams),
  result_hessian(config.hessian_inner.nnz())
  {
  // copy in beta and theta from config
    for(size_t D=0,Dbeta=config.beta_begin;Dbeta<config.beta_end;D++,Dbeta++) {
      parameters[Dbeta] = config.beta[D];
    }
    for(size_t D=0,Dtheta=config.theta_begin;Dtheta<config.theta_end;D++,Dtheta++) {
      parameters[Dtheta] = config.theta[D];
    }

    const CscPattern& hessianIJ = config.hessian_inner;

    const int n   = static_cast<int>(hessianIJ.nrow());
    const int nnz = static_cast<int>(hessianIJ.nnz());

    const int* p = hessianIJ.p.data();
    const int* r = hessianIJ.i.data();

    h_index_upper.resize(nnz);
    std::iota(h_index_upper.begin(), h_index_upper.end(), 0);

    using SpMatI = Eigen::SparseMatrix<int, Eigen::ColMajor, int>;
    Eigen::Map<const SpMatI> Hmap(n, n, nnz, p, r, h_index_upper.data());

  // Materialize full symmetric mapping
    Htemplate = Hmap.selfadjointView<Eigen::Upper>();
    Htemplate.makeCompressed();
  }

public:
// -------------------------
// trustOptim interface bits
// -------------------------

  int get_nvars() const {
    return (int) Htemplate.cols();
  }

  int get_nnz() const {
    return (int) Htemplate.nonZeros();
  }

// f only
template <class DerivedX>
  void get_f(const Eigen::MatrixBase<DerivedX> &x, double &f) {

    for(size_t D=0, Dgamma=config.gamma_begin; Dgamma<config.gamma_end; ++D, ++Dgamma) {
      parameters[Dgamma] = x[D];
    }

    double f_sum = 0.0;
    const int end_group = (int)config.Sgroups.size();

#pragma omp parallel for schedule(dynamic) reduction(+:f_sum)
    for (int kk = 0; kk < end_group; ++kk) {
      const size_t Dgroup = config.Sgroups[(size_t)kk];
      f_sum += fval_api(parameters, Dgroup, adPack);
    }
// return minus log likelihood
    f = -f_sum;
  }

// f and g
template <class DerivedX, class DerivedG>
  void get_fdf(const Eigen::MatrixBase<DerivedX> &x, double &f, Eigen::MatrixBase<DerivedG> &g) {

    for(size_t D=0,Dgamma=config.gamma_begin;Dgamma<config.gamma_end;D++,Dgamma++) {
      parameters[Dgamma] = x[D];
    }

    f = 0.0;
    auto& gout = g.derived();
    gout.setZero();
    const CscPattern& pattern = config.match.grad_inner;

  #pragma omp parallel
    {
      double f_local=0.0;
      CPPAD_TESTVECTOR(double) result_grad_local(config.Ngamma, 0.0);
      CPPAD_TESTVECTOR(double) result_grad_single(config.Ngamma);

      const int end_group = (int) config.Sgroups.size();
    #pragma omp for
      for (int kk = 0; kk < end_group; ++kk) {
        const size_t Dgroup = config.Sgroups[kk];
        const size_t Nhere = pattern.p[Dgroup+1] - pattern.p[Dgroup];
        result_grad_single.resize(Nhere);

        f_local += grad_api(parameters, Dgroup, adPack, true, result_grad_single);
        for(size_t Dlocal=0,Di=pattern.p[Dgroup];Dlocal < Nhere; Dlocal++,Di++) {
          result_grad_local[pattern.i[Di]] += result_grad_single[Dlocal];
        }
      }
    #pragma omp critical 
      {
              // minus so we return minus log likelihood

        f -= f_local;
        for(size_t D=0;D < config.Ngamma; D++) {
          gout[(Eigen::Index) D] -= result_grad_local[D];
        }
      }
  } // parallel

}

// f, g, and H
template <class DerivedX, class DerivedG>
void get_fdfh(const Eigen::MatrixBase<DerivedX> &x,
  double &f,
  Eigen::MatrixBase<DerivedG> &g,
  Eigen::SparseMatrix<double> &H) {

  for(size_t D=0,Dgamma=config.gamma_begin;Dgamma<config.gamma_end;D++,Dgamma++) {
    parameters[Dgamma] = x[D];
  }
  f = 0.0;
  auto& gout = g.derived();
  if (gout.size() > 0) {
    gout.setZero();
  }
  std::fill(result_hessian.begin(), result_hessian.end(), 0.0);

  const CscPattern& pattern_grad = config.match.grad_inner;
  const CscPattern& pattern_hess = config.match.hessian_inner;
  const size_t Nhess = result_hessian.size();

  #pragma omp parallel
  {
    double f_local=0.0;
    CPPAD_TESTVECTOR(double) result_grad_local(config.Ngamma, 0.0);
    CPPAD_TESTVECTOR(double) result_grad_single(config.Ngamma);
    CPPAD_TESTVECTOR(double) result_hess_local(Nhess, 0.0);
    CPPAD_TESTVECTOR(double) result_hess_single(Nhess);

    const int end_group = (int) config.Sgroups.size();
    #pragma omp for schedule(dynamic)
    for (int kk = 0; kk < end_group; ++kk) {
      const size_t Dgroup = config.Sgroups[kk];
      const size_t NhereGrad = pattern_grad.p[Dgroup+1] - pattern_grad.p[Dgroup];
      result_grad_single.resize(NhereGrad);

      const std::size_t start = pattern_hess.p[Dgroup];
      const std::size_t end   = pattern_hess.p[Dgroup + 1];
      const std::size_t Nhere = end - start;

      result_hess_single.resize(Nhere);

      f_local += hess_api(parameters, Dgroup, adPack, true, result_grad_single, result_hess_single);
      for(size_t D=0,Di=pattern_grad.p[Dgroup];D < NhereGrad; D++,Di++) {
        result_grad_local[pattern_grad.i[Di]] += result_grad_single[D];
      }
      for(size_t D=0,Di=start;D < Nhere; D++,Di++) {
        result_hess_local[pattern_hess.i[Di]] += result_hess_single[D];
      }
    }
    #pragma omp critical 
    {
      // minus so we return minus log likelihood
      f -= f_local;
      for(size_t D=0;D < config.Ngamma; D++) {
        gout[D] -= result_grad_local[D];
      }
      for(size_t D=0;D < Nhess; D++) {
        result_hessian[D] -= result_hess_local[D];
      }
    }
  } // parallel


  const Eigen::Index Nout = H.nonZeros();

#ifdef DEBUG
// H must already have the same sparsity structure/order as Htemplate.cast<double>()
  if (Nout != Htemplate.nonZeros()) {
    Rcpp::Rcout << "Hnonzeros " << H.nonZeros() << " Htemplate " << Htemplate.nonZeros() << "\n";
    Rcpp::stop("H structure mismatch: nonZeros differ");
  }
#endif
  
  double* Hx = H.valuePtr();
  const int* map = Htemplate.valuePtr();

  for(Eigen::Index D=0; D<Nout;D++) {
    const size_t inUpper = (size_t) map[D];
    Hx[D] = result_hessian[inUpper];
  }
}

template <class DerivedX>
void get_hessian(
  const Eigen::MatrixBase<DerivedX>& x,
  Eigen::SparseMatrix<double>& H
  ) {
  double f_dummy;
  Eigen::VectorXd g_dummy(0);

// Make sure H has the correct sparsity
  if (H.nonZeros() != Htemplate.nonZeros()) {
    H = Htemplate.cast<double>();
  }

// Reuse your full evaluator
  get_fdfh(x, f_dummy, g_dummy, H);
}

template <class DerivedX, class DerivedH>
void get_hessian(const Eigen::MatrixBase<DerivedX>& x,
                 Eigen::SparseMatrixBase<DerivedH>& Hbase)
{
  DerivedH& H = Hbase.derived();

  // Ensure correct sparsity/order (so your valuePtr() mapping is valid)
  if (H.rows() != Htemplate.rows() ||
      H.cols() != Htemplate.cols() ||
      H.nonZeros() != Htemplate.nonZeros())
  {
    H = Htemplate.cast<double>();   // copies structure + allocates values
    H.makeCompressed();
  }

  double f_dummy = 0.0;
  Eigen::VectorXd g_dummy(0);       // sized 0: ok if get_fdfh doesn't write g blindly

  get_fdfh(x, f_dummy, g_dummy, H);
}

template <class DerivedX, class DerivedG>
void get_df(const Eigen::MatrixBase<DerivedX>& x,
  Eigen::MatrixBase<DerivedG>& g)
{
  double f_dummy;
  get_fdf(x, f_dummy, g);
}
};


#endif
