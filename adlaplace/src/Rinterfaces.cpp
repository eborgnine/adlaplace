#include <Rcpp.h>
#include <cppad/cppad.hpp>

#include <Rinternals.h>
#include "adlaplace/creators/R_interfaces.hpp"

double jointLogDens(const Rcpp::NumericVector& x, SEXP backendContext, SEXP Sgroups);
Rcpp::NumericVector grad(const Rcpp::NumericVector& x, SEXP backendContext, bool inner, SEXP Sgroups);
Rcpp::S4 hess(const Rcpp::NumericVector& x, SEXP backendContext, bool inner, SEXP Sgroups);
Rcpp::NumericVector traceHinvT(
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  SEXP backendContext,
  SEXP Sgroups
);

extern "C" SEXP adlaplace_grad_c(SEXP xSEXP, SEXP backendContextSEXP, SEXP innerSEXP, SEXP SgroupsSEXP) {
  BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  const Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xSEXP);
  const bool inner = Rcpp::as<bool>(innerSEXP);
  const Rcpp::NumericVector out = grad(x, backendContextSEXP, inner, SgroupsSEXP);
  return Rcpp::wrap(out);
  END_RCPP
}

extern "C" SEXP adlaplace_hess_c(SEXP xSEXP, SEXP backendContextSEXP, SEXP innerSEXP, SEXP SgroupsSEXP) {
  BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  const Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xSEXP);
  const bool inner = Rcpp::as<bool>(innerSEXP);
  const Rcpp::S4 out = hess(x, backendContextSEXP, inner, SgroupsSEXP);
  return Rcpp::wrap(out);
  END_RCPP
}

extern "C" SEXP adlaplace_jointLogDens_c(SEXP xSEXP, SEXP backendContextSEXP, SEXP SgroupsSEXP) {
  BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  const Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xSEXP);
  const double out = jointLogDens(x, backendContextSEXP, SgroupsSEXP);
  return Rcpp::wrap(out);
  END_RCPP
}

extern "C" SEXP adlaplace_traceHinvT_c(
  SEXP xSEXP,
  SEXP backendContextSEXP,
  SEXP LinvPtSEXP,
  SEXP LinvPtColumnsSEXP,
  SEXP SgroupsSEXP
) {
  BEGIN_RCPP
  Rcpp::RNGScope rcpp_rngScope_gen;
  const Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(xSEXP);
  const Rcpp::S4 LinvPt(LinvPtSEXP);
  const Rcpp::S4 LinvPtColumns(LinvPtColumnsSEXP);
  const Rcpp::NumericVector out = traceHinvT(
    x,
    LinvPt,
    LinvPtColumns,
    backendContextSEXP,
    SgroupsSEXP
  );
  return Rcpp::wrap(out);
  END_RCPP
}

//' C++ backend entry points
//'
//' Low-level C++ entry points exposed to R via Rcpp.
//' These create and operate on an opaque backend handle (external pointer)
//' used to evaluate objective, gradient, and Hessian values.
//'
//' @param data An R list containing model data objects required by the backend
//'   (used by \code{getAdFun()}).
//' @param config An R list of configuration options required by the backend
//'   (used by \code{getAdFun()}).
//' @param x Numeric parameter vector of length \code{Nparams}.
//' @param backendContext External pointer returned by \code{getAdFun()}.
//' @param inner Logical scalar. If \code{TRUE}, evaluate inner-\eqn{\gamma}
//'   derivatives; otherwise evaluate outer/full derivatives.
//' @param Sgroups Optional integer vector of 0-based group indices to evaluate.
//'   If omitted, uses all groups \code{0:(Ngroups-1)}.
//' @param LinvPt Sparse \code{dgCMatrix} for columns of
//'   \eqn{P^\top L^{-1} D^{-1/2}} (or equivalent) used in trace contractions.
//' @param LinvPtColumns Sparse \code{ngCMatrix}/\code{dgCMatrix} mapping
//'   selected columns of \code{LinvPt} to each group.
//'
//' @return
//' \itemize{
//'   \item \code{getAdFun}: external pointer handle with backend state.
//'   \item \code{jointLogDens}: scalar objective value summed over groups.
//'   \item \code{grad}: numeric gradient vector.
//'   \item \code{hess}: sparse symmetric Hessian as a Matrix
//'     \code{dsCMatrix} object.
//'   \item \code{traceHinvT}: numeric vector of third-derivative contractions.
//' }
//'
//' @details
//' The external pointer returned by \code{getAdFun()} is opaque and not
//' user-modifiable. It may hold substantial memory (AD tapes, sparsity maps,
//' work caches). Do not save it across R sessions.
//'
//' @name adlaplace_cpp

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
SEXP getAdFun(
  Rcpp::List data,
  Rcpp::List config)
{
  return getAdFun_h(data, config);
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
double jointLogDens(
  const Rcpp::NumericVector& x,
  SEXP backendContext,
  SEXP Sgroups = R_NilValue) {
  adlaplace_adpack_handle* h = get_handle(backendContext);

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  double total = 0.0;
  for (size_t g : groups) {
    double fg = 0.0;
    int gi = static_cast<int>(g);
    const int rc = h->api->f(h->ctx, &gi, x.begin(), &fg);
    if (rc != 0) {
      Rcpp::stop("backend api->f failed for group %d with code %d", gi, rc);
    }
    total += fg;
  }
  return total;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
  const Rcpp::NumericVector& x,
  SEXP backendContext,
  const bool inner = false,
  SEXP Sgroups = R_NilValue) {
  adlaplace_adpack_handle* h = get_handle(backendContext);
  if (!h->api->f_grad) {
    Rcpp::stop("backendContext api->f_grad is NULL");
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  Rcpp::NumericVector grad_out(Nparams, 0.0);
  double f_total = 0.0;
  const bool inner_local = inner;
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  for (size_t g : groups) {
    int gi = static_cast<int>(g);
    const int rc = h->api->f_grad(
      h->ctx, &gi, x.begin(), &inner_local, &f_total, grad_out.begin()
    );
    if (rc != 0) {
      Rcpp::stop("backend api->f_grad failed for group %d with code %d", gi, rc);
    }
  }
  return grad_out;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::S4 hess(
  const Rcpp::NumericVector& x,
  SEXP backendContext,
  const bool inner = false,
  SEXP Sgroups = R_NilValue) {
  adlaplace_adpack_handle* h = get_handle(backendContext);
  if (!h->api->f_grad_hess) {
    Rcpp::stop("backendContext api->f_grad_hess is NULL");
  }
  if (!h->api->get_hessian) {
    Rcpp::stop("backendContext api->get_hessian is NULL");
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  const bool inner_local = inner;
  const int* hess_p = NULL;
  const int* hess_i = NULL;
  size_t hess_p_len = 0;
  size_t hess_i_len = 0;
  const int rc_hessian = h->api->get_hessian(
    h->ctx, &inner_local, &hess_p, &hess_p_len, &hess_i, &hess_i_len
  );
  if (rc_hessian != 0) {
    Rcpp::stop("backend api->get_hessian failed with code %d", rc_hessian);
  }

  Rcpp::NumericVector hess_out(hess_i_len, 0.0);
  Rcpp::NumericVector grad_scratch(Nparams, 0.0);
  double f_total = 0.0;
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  for (size_t g : groups) {
    int gi = static_cast<int>(g);
    const int rc = h->api->f_grad_hess(
      h->ctx,
      &gi,
      x.begin(),
      &inner_local,
      &f_total,
      grad_scratch.begin(),
      hess_out.begin()
    );
    if (rc != 0) {
      Rcpp::stop("backend api->f_grad_hess failed for group %d with code %d", gi, rc);
    }
  }

  const int ncol = static_cast<int>(hess_p_len > 0 ? hess_p_len - 1 : 0);
  Rcpp::IntegerVector p_out(hess_p_len);
  Rcpp::IntegerVector i_out(hess_i_len);
  for (size_t k = 0; k < hess_p_len; ++k) p_out[k] = static_cast<int>(hess_p[k]);
  for (size_t k = 0; k < hess_i_len; ++k) i_out[k] = static_cast<int>(hess_i[k]);

  Rcpp::Environment matrix_ns = Rcpp::Environment::namespace_env("Matrix");
  (void)matrix_ns;
  Rcpp::S4 out("dsCMatrix");
  out.slot("i") = i_out;
  out.slot("p") = p_out;
  out.slot("x") = hess_out;
  out.slot("Dim") = Rcpp::IntegerVector::create(ncol, ncol);
  out.slot("uplo") = Rcpp::String("U");
  return out;
}

//' @rdname adlaplace_cpp
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
  const Rcpp::NumericVector& x,
  const Rcpp::S4& LinvPt,
  const Rcpp::S4& LinvPtColumns,
  SEXP backendContext,
  SEXP Sgroups = R_NilValue
) {
  adlaplace_adpack_handle* h = get_handle(backendContext);
  if (!h->api->trace_hinv_t) {
    Rcpp::stop("backendContext api->trace_hinv_t is NULL");
  }

  size_t Nparams = 0, Ngroups = 0, Nbeta = 0, Ngamma = 0, Ntheta = 0;
  const int rc_sizes = h->api->get_sizes(
    h->ctx, &Nparams, &Ngroups, &Nbeta, &Ngamma, &Ntheta
  );
  if (rc_sizes != 0) {
    Rcpp::stop("backend api->get_sizes failed with code %d", rc_sizes);
  }
  if (static_cast<size_t>(x.size()) != Nparams) {
    Rcpp::stop("x has length %d but expected Nparams=%d", x.size(), (int)Nparams);
  }

  Rcpp::IntegerVector LinvPt_p = LinvPt.slot("p");
  Rcpp::IntegerVector LinvPt_i = LinvPt.slot("i");
  Rcpp::NumericVector LinvPt_x = LinvPt.slot("x");
  Rcpp::IntegerVector LinvPt_Dim = LinvPt.slot("Dim");

  Rcpp::IntegerVector LinvPtColumns_p = LinvPtColumns.slot("p");
  Rcpp::IntegerVector LinvPtColumns_i = LinvPtColumns.slot("i");

  const size_t LinvPt_ncol = static_cast<size_t>(LinvPt_Dim[1]);
  const size_t LinvPt_p_len = static_cast<size_t>(LinvPt_p.size());
  const size_t LinvPt_i_len = static_cast<size_t>(LinvPt_i.size());
  const size_t LinvPt_x_len = static_cast<size_t>(LinvPt_x.size());
  const size_t LinvPtColumns_p_len = static_cast<size_t>(LinvPtColumns_p.size());
  const size_t LinvPtColumns_i_len = static_cast<size_t>(LinvPtColumns_i.size());

  Rcpp::NumericVector trace_out(Nparams, 0.0);
  const Rcpp::IntegerVector Sgroups_vec = (Sgroups == R_NilValue)
    ? Rcpp::IntegerVector()
    : Rcpp::as<Rcpp::IntegerVector>(Sgroups);
  const std::vector<size_t> groups = resolve_groups(Ngroups, Sgroups_vec);

  for (size_t g : groups) {
    int gi = static_cast<int>(g);
    const int rc = h->api->trace_hinv_t(
      h->ctx,
      &gi,
      x.begin(),
      LinvPt_p.begin(),
      LinvPt_i.begin(),
      LinvPt_x.begin(),
      LinvPt_ncol,
      LinvPt_p_len,
      LinvPt_i_len,
      LinvPt_x_len,
      LinvPtColumns_p.begin(),
      LinvPtColumns_i.begin(),
      LinvPtColumns_p_len,
      LinvPtColumns_i_len,
      trace_out.begin()
    );
    if (rc != 0) {
      Rcpp::stop("backend api->trace_hinv_t failed for group %d with code %d", gi, rc);
    }
  }

  return trace_out;
}
