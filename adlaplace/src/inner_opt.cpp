
//' Inner optimization over gamma using trust-region CG (sparse)
//'
//' Runs the inner optimization problem (typically over \eqn{\gamma}) using the
//' trustOptim sparse trust-region Conjugate Gradient solver. This function
//' evaluates the objective, gradient, and Hessian through the pre-built AD pack
//' (external pointer) and returns the solution along with curvature information.
//'
//' @param x Numeric full parameter vector of length \code{Nparams}.
//' @param parameters Numeric vector of fixed outer parameters
//'   (\code{beta}, \code{theta}; length \code{Nbeta+Ntheta}) used by
//'   \code{inner_opt()}.
//' @param gamma Numeric vector of starting values for inner parameters
//'   (\code{gamma}; length \code{Ngamma}) used by \code{inner_opt()}.
//' @param ad_fun \code{ad_fun} S4 object from \code{ad_fun(ad_fun_ptr)}.
//' @param verbose Logical; if \code{TRUE}, print threads, shards, and parameter sizes.
//' @param control List of trust-region control parameters for
//'   \code{inner_opt()} (see \pkg{trustOptim}).
//' @param deriv Logical: if \code{TRUE}, return full outer gradient and Hessian at the
//'   inner solution; if \code{FALSE}, return inner quantities only (default).
//'
//' @return
//'   \item{\code{inner_opt()}}{Returns \code{log_lik}, \code{neg_log_lik}
//'   (Laplace profile likelihood and its negation), \code{parameters} (outer
//'   \code{beta} and \code{theta} passed in), \code{full_parameters} (outer plus
//'   inner \code{gamma} at the mode), \code{fval} (inner objective:
//'   negative log joint density at \eqn{\hat\gamma}), \code{solution},
//'   \code{gradient} (list \code{inner}, \code{outer}; \code{outer} empty when
//'   \code{deriv=FALSE}), \code{hessian} (list \code{inner}, \code{outer},
//'   \code{chol_inner}, \code{half_log_det}; when \code{deriv=TRUE} also
//'   \code{half_H_inv}, \code{H_inv}, and \code{trace3}),
//'   \code{iterations}, \code{status}, \code{trust.radius}, \code{method}.
//'   Objective and derivatives use the **negative log-density** convention.}
//'
//' @details
//' This calls the sparse method from the \code{TrustOptim} package via the Cpp interface.
//' \code{inner_opt()} negates tape log-density values in C++ for minimization.
//'
//' @name innerOpt

// Standard C
#include <cmath>
#include <cstdio>
#include <cstdlib>
// Eigen
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

// Rcpp
#include <Rcpp.h>

// trustOptim
#include <CG-sparse.h>

// Local
#include "adlaplace/math/constants.hpp"
#include "adlaplace/creators/rviews.hpp"
#include "adlaplace/ompad.hpp"
#include "adlaplace/runtime/thread_groups.hpp"
#include "adlaplace/runtime/thread_affinity_debug.hpp"
#include "chol_update.hpp"
#include "trace_hinv_t_runtime.hpp"
#include "trustOptimWrappers.hpp"
#include "trustOptimControl.hpp"

// Convert Eigen sparse matrix to Matrix::dgCMatrix S4 object directly
Rcpp::S4 eigen_to_dgCMatrix(const Eigen::SparseMatrix<double>& M) {
	const Eigen::Index nrow = M.rows();
	const Eigen::Index ncol = M.cols();
	if (nrow == 0 || ncol == 0) {
		Rcpp::S4 mat("dgCMatrix");
		mat.slot("p") = Rcpp::IntegerVector::create(0);
		mat.slot("i") = Rcpp::IntegerVector();
		mat.slot("x") = Rcpp::NumericVector();
		mat.slot("Dim") = Rcpp::IntegerVector::create(0, 0);
		return mat;
	}
	const Eigen::Index nnz = M.nonZeros();

	Rcpp::IntegerVector i(nnz);
	Rcpp::IntegerVector p(ncol + 1);
	Rcpp::NumericVector x(nnz);

	// Copy data from Eigen (already in CSC format, 0-based)
	std::copy(M.innerIndexPtr(), M.innerIndexPtr() + nnz, i.begin());
	std::copy(M.outerIndexPtr(), M.outerIndexPtr() + ncol + 1, p.begin());
	std::copy(M.valuePtr(), M.valuePtr() + nnz, x.begin());

	Rcpp::S4 mat("dgCMatrix");
	mat.slot("i") = i;
	mat.slot("p") = p;
	mat.slot("x") = x;
	mat.slot("Dim") = Rcpp::IntegerVector::create(static_cast<int>(nrow), static_cast<int>(ncol));
	return mat;
}

// Unit-lower triangular factor L from chol_inner CSC pattern (0-based p, i).
Rcpp::S4 csc_to_dtCMatrix_lower(
	const std::vector<int>& p,
	const std::vector<int>& i,
	const std::vector<double>& x,
	int n)
{
	if (static_cast<size_t>(n + 1) != p.size()) {
		Rcpp::stop("csc_to_dtCMatrix_lower: p length must be nrow + 1");
	}
	if (i.size() != x.size()) {
		Rcpp::stop("csc_to_dtCMatrix_lower: i and x length must match");
	}

	Rcpp::S4 mat("dtCMatrix");
	mat.slot("p") = Rcpp::IntegerVector(p.begin(), p.end());
	mat.slot("i") = Rcpp::IntegerVector(i.begin(), i.end());
	mat.slot("x") = Rcpp::NumericVector(x.begin(), x.end());
	mat.slot("Dim") = Rcpp::IntegerVector::create(n, n);
	mat.slot("uplo") = Rcpp::String("L");
	mat.slot("diag") = Rcpp::String("N");
	return mat;
}

Rcpp::S4 csc_to_dgCMatrix(
	const std::vector<int>& p,
	const std::vector<int>& i,
	const std::vector<double>& x,
	int n)
{
	if (static_cast<size_t>(n + 1) != p.size()) {
		Rcpp::stop("csc_to_dgCMatrix: p length must be nrow + 1");
	}
	if (i.size() != x.size()) {
		Rcpp::stop("csc_to_dgCMatrix: i and x length must match");
	}

	Rcpp::S4 mat("dgCMatrix");
	mat.slot("p") = Rcpp::IntegerVector(p.begin(), p.end());
	mat.slot("i") = Rcpp::IntegerVector(i.begin(), i.end());
	mat.slot("x") = Rcpp::NumericVector(x.begin(), x.end());
	mat.slot("Dim") = Rcpp::IntegerVector::create(n, n);
	return mat;
}

// Numeric LDL from chol_update (x on L1 pattern, d on diagonal); pattern from chol_inner_list.
Rcpp::List chol_inner_numeric(
	const ad_fun& backend,
	const std::vector<double>& x,
	const std::vector<double>& d,
	const std::vector<double>& linv_x,
	bool include_linv)
{
	const CholPattern& pat = backend.chol_pattern;

	Rcpp::List out = Rcpp::List::create(
		Rcpp::Named("L1") = csc_to_dtCMatrix_lower(pat.L1_p, pat.L1_i, x, pat.n),
		Rcpp::Named("D") = Rcpp::NumericVector(d.begin(), d.end()),
		Rcpp::Named("perm") = Rcpp::IntegerVector(pat.perm.begin(), pat.perm.end()),
		Rcpp::Named("perm_inv") = Rcpp::IntegerVector(pat.perm_inv.begin(), pat.perm_inv.end())
	);

	if (include_linv) {
		out["Linv"] = csc_to_dgCMatrix(
			pat.Linv_p,
			pat.Linv_i,
			linv_x,
			pat.n
		);
	}

	return out;
}

struct InnerOptResult {
	double fval = NA_REAL;
	double log_lik = NA_REAL;
	double neg_log_lik = NA_REAL;
	Eigen::VectorXd solution;
	Eigen::VectorXd parameters;
	Eigen::VectorXd full_parameters;
	Eigen::VectorXd grad_inner;
	Eigen::VectorXd grad_outer;
	Eigen::SparseMatrix<double> hessian_inner;
	Eigen::SparseMatrix<double> hessian_outer;
	std::vector<double> x_out;
	std::vector<double> d_out;
	std::vector<double> linv_x_out;
	std::vector<double> half_h_inv_x_out;
	std::vector<double> h_inv_x_out;
	std::vector<double> trace3;
	int iterations = NA_INTEGER;
	MB_Status status = SUCCESS;
	double trust_radius = NA_REAL;
	bool chol_ok = false;
	double half_log_det = NA_REAL;
};

InnerOptResult inner_opt(
	const std::vector<double>& parameters,
	const std::vector<double>& gamma,
	ad_fun& backend,
	const TrustControl& control,
	bool deriv,
	bool verbose)
{
	adlaplace_require_owner_threads_assigned(backend);
	const std::vector<std::vector<size_t>> thread_groups = thread_groups_from_backend(backend);
	const int num_threads = static_cast<int>(thread_groups.size());
	using Tvec = Eigen::VectorXd;
	using THess = Eigen::SparseMatrix<double>;
	using TPreLLt = Eigen::SimplicialLLT<THess>;

	const std::size_t n_beta = static_cast<std::size_t>(backend.sizes.named("beta"));
	const std::size_t n_gamma = static_cast<std::size_t>(backend.sizes.named("gamma"));
	const std::size_t n_theta = static_cast<std::size_t>(backend.sizes.named("theta"));
	const std::size_t n_params = n_beta + n_gamma + n_theta;
	const std::size_t gamma_begin = n_beta;
	const std::size_t theta_begin = n_beta + n_gamma;

	if (parameters.size() != n_beta + n_theta) {
		Rcpp::stop(
			"parameters has length %d but expected Nbeta+Ntheta=%d",
			static_cast<int>(parameters.size()),
			static_cast<int>(n_beta + n_theta)
		);
	}
	if (gamma.size() != n_gamma) {
		Rcpp::stop(
			"gamma has length %d but expected Ngamma=%d",
			static_cast<int>(gamma.size()),
			static_cast<int>(n_gamma)
		);
	}
	if (n_gamma == 0) {
		Rcpp::stop("Ngamma must be > 0");
	}

	if (verbose) {
		Rcpp::Rcout << "inner_opt: threads = " << num_threads
		            << ", shards = " << backend.fun.size()
		            << ", params = " << n_params
		            << " (beta = " << n_beta
		            << ", gamma = " << n_gamma
		            << ", theta = " << n_theta << ")\n";
	}
	if (verbose && adlaplace_debug_enabled()) {
		Rcpp::Rcout << "inner_opt: DEBUG build (thread-affinity checks active)\n";
	}

	Tvec gamma_start(static_cast<Eigen::Index>(n_gamma));
	Tvec solution(static_cast<Eigen::Index>(n_gamma));
	Tvec fullParams(static_cast<Eigen::Index>(n_params));
	Tvec grad(static_cast<Eigen::Index>(n_gamma));

	std::copy(gamma.begin(), gamma.end(), gamma_start.data());

	std::vector<double> params_init(n_params);
	for (std::size_t d = 0; d < n_beta; ++d) {
		params_init[d] = parameters[d];
		fullParams[static_cast<Eigen::Index>(d)] = parameters[d];
	}
	for (std::size_t d = 0; d < n_gamma; ++d) {
		params_init[gamma_begin + d] = gamma[d];
	}
	for (std::size_t d = 0; d < n_theta; ++d) {
		params_init[theta_begin + d] = parameters[n_beta + d];
		fullParams[static_cast<Eigen::Index>(theta_begin + d)] = parameters[n_beta + d];
	}

	AD_Func_Opt funObj(backend, params_init, true, num_threads, &thread_groups);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();

	Eigen::SparseMatrix<double> Houter;
	Tvec gradOuter;

	double fval = NA_REAL;
	double radius = NA_REAL;
	int iterations = NA_INTEGER;
	MB_Status status = SUCCESS;

	InnerOptResult out;
	out.x_out.assign(backend.chol_pattern.L1_i.size(), 0.0);
	out.d_out.assign(H.rows(), 0.0);
	if (deriv) {
		out.linv_x_out.assign(backend.chol_pattern.Linv_i.size(), 0.0);
	}

	double log_det = NA_REAL;

	{
		CppadParallelScope parallel_scope(static_cast<std::size_t>(num_threads));

		Trust_CG_Sparse<Tvec, AD_Func_Opt, THess, TPreLLt> opt(
			funObj, gamma_start,
			control.rad, control.min_rad, control.tol, control.prec,
			control.report_freq, control.report_level, control.header_freq,
			control.report_precision,
			control.maxit, control.contract_factor, control.expand_factor,
			control.contract_threshold, control.expand_threshold_rad,
			control.expand_threshold_ap, control.function_scale_factor,
			control.precond_refresh_freq, control.precond_ID, control.trust_iter
		);

		opt.run();
		status = opt.get_current_state(solution, fval, grad, H, iterations, radius);

		const double grad_l2_sq = grad.squaredNorm();
		if (grad_l2_sq > 10.0) {
			if (verbose && adlaplace_debug_enabled()) {
				Rcpp::Rcout << "restarting with shrunk gamma\n";
			}

			gamma_start = gamma_start.cwiseMax(-0.1).cwiseMin(0.1);
			Trust_CG_Sparse<Tvec, AD_Func_Opt, THess, TPreLLt> opt_retry(
				funObj, gamma_start,
				control.rad, control.min_rad, control.tol, control.prec,
				control.report_freq, control.report_level, control.header_freq,
				control.report_precision,
				control.maxit, control.contract_factor, control.expand_factor,
				control.contract_threshold, control.expand_threshold_rad,
				control.expand_threshold_ap, control.function_scale_factor,
				control.precond_refresh_freq, control.precond_ID, control.trust_iter
			);
			opt_retry.run();
			status = opt_retry.get_current_state(solution, fval, grad, H, iterations, radius);
		}

		H.makeCompressed();
		log_det = chol_update(
			H,
			backend.chol_pattern.perm,
			backend.chol_pattern.L1_p,
			backend.chol_pattern.L1_i,
			out.x_out,
			out.d_out
		);
		
		for (std::size_t d = 0; d < n_gamma; ++d) {
			fullParams[static_cast<Eigen::Index>(gamma_begin + d)] = solution[static_cast<Eigen::Index>(d)];
		}

		if (deriv) {
			// Outer grad/Hessian in the same CppAD parallel session as inner TR.
			if (verbose && adlaplace_debug_enabled()) {
				Rcpp::Rcout << "inner_opt: trust region done; outer get_fdfh next ("
				            << num_threads << " threads)\n";
			}
			AD_Func_Opt funObjOuter(
				backend, params_init, false, num_threads, &thread_groups);

			gradOuter = Tvec(static_cast<Eigen::Index>(n_params));
			Houter = funObjOuter.Htemplate.cast<double>();
			funObjOuter.get_fdfh(fullParams, fval, gradOuter, Houter);

			if (verbose && adlaplace_debug_enabled()) {
				Rcpp::Rcout << "inner_opt: outer get_fdfh done\n";
			}
			adlaplace_debug_raise_if_any("inner_opt outer get_fdfh");

				linv_update(
					backend.chol_pattern.L1_p,
					backend.chol_pattern.L1_i,
					out.x_out,
					backend.chol_pattern.Linv_p,
					backend.chol_pattern.Linv_i,
					out.linv_x_out
				);

				half_h_inv_update(
						backend.chol_pattern.Linv_p,
						backend.chol_pattern.Linv_i,
						out.linv_x_out,
						out.d_out,
						backend.chol_pattern.perm_inv,
						backend.chol_pattern.half_H_inv_p,
						backend.chol_pattern.half_H_inv_i,
						out.half_h_inv_x_out
					);

					h_inv_update(
						backend.chol_pattern.half_H_inv_p,
						backend.chol_pattern.half_H_inv_i,
						out.half_h_inv_x_out,
						backend.chol_pattern.H_inv_p,
						backend.chol_pattern.H_inv_i,
						out.h_inv_x_out
					);

					if (verbose && adlaplace_debug_enabled()) {
						Rcpp::Rcout << "inner_opt: trace_hinv_t next ("
						            << num_threads << " threads)\n";
					}

					const std::vector<double> x_vec(
						fullParams.data(),
						fullParams.data() + static_cast<Eigen::Index>(n_params)
					);
					out.trace3 = adlaplace_trace::trace_hinv_t_parallel(
						backend,
						x_vec,
						backend.chol_pattern.half_H_inv_p,
						backend.chol_pattern.half_H_inv_i,
						out.half_h_inv_x_out,
						n_gamma,
						backend.chol_pattern.trace_columns_p,
						backend.chol_pattern.trace_columns_i,
						verbose
					);

				}
	}

	if (verbose && adlaplace_debug_enabled()) {
		Rcpp::Rcout << "inner_opt: parallel block ended\n";
	}
	adlaplace_debug_raise_if_any("inner_opt after parallel block");

	out.fval = fval;
	out.solution = solution;
	out.parameters = Eigen::VectorXd(static_cast<Eigen::Index>(n_beta + n_theta));
	for (std::size_t d = 0; d < n_beta; ++d) {
		out.parameters[static_cast<Eigen::Index>(d)] = parameters[d];
	}
	for (std::size_t d = 0; d < n_theta; ++d) {
		out.parameters[static_cast<Eigen::Index>(n_beta + d)] = parameters[n_beta + d];
	}
	out.full_parameters = fullParams;
	out.iterations = iterations;
	out.status = status;
	out.trust_radius = radius;

	if (R_finite(log_det)) {
		out.half_log_det = 0.5 * log_det;
		out.log_lik = -out.fval + static_cast<double>(n_gamma) * ONEHALFLOGTWOPI - out.half_log_det;
		out.neg_log_lik = -out.log_lik;
		out.chol_ok = true;
	} else {
		Rcpp::warning(
			"inner_opt: non-positive or non-finite D; half_log_det not computed"
		);
	}

	out.grad_inner = grad;
	out.hessian_inner = H;
	if (deriv) {
		out.grad_outer = gradOuter;
		out.hessian_outer = Houter;
	} else {
		out.grad_outer.resize(0);
		out.hessian_outer.resize(0, 0);
	}

	if (verbose) {
		Rcpp::Rcout << "inner_opt: finished\n";
	}

	return out;
}


//' @rdname innerOpt
//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
	const Rcpp::NumericVector parameters,
	const Rcpp::NumericVector gamma,
	const Rcpp::S4& ad_fun,
	SEXP control = R_NilValue,
	bool deriv = false,
	bool verbose = false)
{
	try {
		const Rcpp::List control_list(
			Rf_isNull(control) ? Rcpp::List() : Rcpp::as<Rcpp::List>(control)
		);
		const TrustControl control_c(control_list);
		::ad_fun* backend = resolve_ad_fun_laplace(ad_fun);

		std::vector<double> parameters_vec(parameters.begin(), parameters.end());
		std::vector<double> gamma_vec(gamma.begin(), gamma.end());

		const InnerOptResult result = inner_opt(
			parameters_vec,
			gamma_vec,
			*backend,
			control_c,
			deriv,
			verbose
		);

		const Rcpp::List chol_inner = chol_inner_numeric(
			*backend,
			result.x_out,
			result.d_out,
			result.linv_x_out,
			deriv
		);

		Rcpp::List hessian_out = Rcpp::List::create(
			Rcpp::Named("inner") = eigen_to_dgCMatrix(result.hessian_inner),
			Rcpp::Named("outer") = eigen_to_dgCMatrix(result.hessian_outer),
			Rcpp::Named("chol_inner") = chol_inner,
			Rcpp::Named("half_log_det") = Rcpp::wrap(result.half_log_det)
		);
		if (deriv && result.chol_ok && !result.half_h_inv_x_out.empty()) {
			const CholPattern& pat = backend->chol_pattern;
			hessian_out["half_H_inv"] = csc_to_dgCMatrix(
				pat.half_H_inv_p,
				pat.half_H_inv_i,
				result.half_h_inv_x_out,
				pat.n
			);
			hessian_out["H_inv"] = csc_to_dgCMatrix(
				pat.H_inv_p,
				pat.H_inv_i,
				result.h_inv_x_out,
				pat.n
			);
		}
		if (deriv && result.chol_ok && !result.trace3.empty()) {
			hessian_out["trace3"] = Rcpp::wrap(result.trace3);
		}

		return Rcpp::List::create(
			Rcpp::Named("log_lik") = Rcpp::wrap(result.log_lik),
			Rcpp::Named("neg_log_lik") = Rcpp::wrap(result.neg_log_lik),
			Rcpp::Named("parameters") = Rcpp::wrap(result.parameters),
			Rcpp::Named("full_parameters") = Rcpp::wrap(result.full_parameters),
			Rcpp::Named("opt") = Rcpp::List::create(
				Rcpp::Named("fval") = Rcpp::wrap(result.fval),
				Rcpp::Named("solution") = Rcpp::wrap(result.solution),
				Rcpp::Named("iterations") = Rcpp::wrap(result.iterations),
				Rcpp::Named("status") = Rcpp::wrap(std::string(MB_strerror(result.status))),
				Rcpp::Named("trust.radius") = Rcpp::wrap(result.trust_radius),
				Rcpp::Named("method") = Rcpp::wrap("Sparse")
			),
			Rcpp::Named("gradient") = Rcpp::List::create(
				Rcpp::Named("inner") = Rcpp::wrap(result.grad_inner),
				Rcpp::Named("outer") = Rcpp::wrap(result.grad_outer)
			),
			Rcpp::Named("hessian") = hessian_out
		);
	}
	catch (const Rcpp::exception& e) {
		Rcpp::stop("inner_opt failed: %s", e.what());
	}
	catch (...) {
		Rcpp::stop("inner_opt failed with unknown error");
	}
}
