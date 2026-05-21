
//' Inner optimization over gamma using trust-region CG (sparse)
//'
//' Runs the inner optimization problem (typically over \eqn{\gamma}) using the
//' trustOptim sparse trust-region Conjugate Gradient solver. This function
//' evaluates the objective, gradient, and Hessian through the pre-built AD pack
//' (external pointer) and returns the solution along with curvature information.
//'
//' @param x Numeric full parameter vector of length \code{Nparams} used by
//'   \code{all_derivs()}.
//' @param parameters Numeric vector of fixed outer parameters
//'   (\code{beta}, \code{theta}; length \code{Nbeta+Ntheta}) used by
//'   \code{inner_opt()}.
//' @param gamma Numeric vector of starting values for inner parameters
//'   (\code{gamma}; length \code{Ngamma}) used by \code{inner_opt()}.
//' @param ad_fun List returned by \code{getAdFun()} (must contain \code{ad_fun}).
//' @param config Configuration list with model dimensions, groups, and
//'   sparsity information.
//' @param control List of trust-region control parameters for
//'   \code{inner_opt()} (see \pkg{trustOptim}).
//' @param deriv Logical: if \code{TRUE}, return full outer gradient and Hessian at the
//'   inner solution; if \code{FALSE}, return inner quantities only. When missing,
//'   uses \code{config$deriv} or legacy \code{config$inner_only}.
//'
//' @return
//'   \item{\code{all_derivs()}}{List with \code{fval}, \code{gradient}, and
//'   \code{hessian} of the **negative log density** \eqn{-\ell(x)} at full
//'   outer parameters \code{x} (minimization / \pkg{trustOptim} sign). Compare
//'   \code{\link{jointLogDens}}, \code{\link{grad}}, and \code{\link{hessian}}
//'   for \eqn{\ell(x)} and its derivatives at the same point.}
//'   \item{\code{inner_opt()}}{Returns \code{fval}, \code{solution},
//'   \code{gradient} (inner gradient if \code{deriv=FALSE}, otherwise full gradient),
//'   \code{hessian} (inner Hessian if \code{deriv=FALSE}, otherwise full Hessian),
//'   \code{iterations}, \code{status}, \code{trust.radius}, \code{method},
//'   \code{chol_inner} (list \code{L1}, \code{D}, \code{perm}, as \code{Matrix::expand2}).
//'   Objective and derivatives use the same **negative log-density** convention as
//'   \code{all_derivs()}.}
//'
//' @details
//' This calls the sparse method from the \code{TrustOptim} package via the Cpp interface.
//' \code{all_derivs()} and \code{inner_opt()} negate the tape log-density values in
//' C++; \code{\link{jointLogDens}} / \code{\link{grad}} / \code{\link{hessian}} do not.  
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
#include "adlaplace/runtime/rviews.hpp"
#include "adlaplace/ompad.hpp"
#include "chol_update.hpp"
#include "trustOptimWrappers.hpp"
#include "trustOptimControl.hpp"

// Convert Eigen sparse matrix to Matrix::dgCMatrix S4 object directly
Rcpp::S4 eigen_to_dgCMatrix(const Eigen::SparseMatrix<double>& M) {
	const Eigen::Index nrow = M.rows();
	const Eigen::Index ncol = M.cols();
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

// L1 matches Matrix::expand2(); D is numeric diagonal (expand2 uses D@x on ddiMatrix).
Rcpp::List chol_inner_list(
	const CholPattern& chol_pat,
	const std::vector<double>& x_out,
	const std::vector<double>& d_out)
{
	return Rcpp::List::create(
		Rcpp::Named("L1") = csc_to_dtCMatrix_lower(chol_pat.p, chol_pat.i, x_out, chol_pat.n),
		Rcpp::Named("D") = Rcpp::NumericVector(d_out.begin(), d_out.end()),
		Rcpp::Named("perm") = Rcpp::IntegerVector(chol_pat.perm.begin(), chol_pat.perm.end())
	);
}

struct InnerOptResult {
	double fval = NA_REAL;
	double log_lik = NA_REAL;
	Eigen::VectorXd solution;
	Eigen::VectorXd full_parameters;
	Eigen::VectorXd gradient;
	Eigen::SparseMatrix<double> hessian;
	std::vector<double> x_out;
	std::vector<double> d_out;
	int iterations = NA_INTEGER;
	MB_Status status = SUCCESS;
	double trust_radius = NA_REAL;
	bool chol_ok = false;
	double half_log_det = NA_REAL;
};

InnerOptResult all_derivs(
	const std::vector<double>& x,
	ad_groups& groups,
	const Config& config)
{
	const int num_threads = config.num_threads > 0 ? config.num_threads : 1;

	AD_Func_Opt funObj(
		groups,
		x,
		false,           // inner=false
		num_threads
	);

	const Eigen::Index nvars = static_cast<Eigen::Index>(funObj.get_nvars());
	Eigen::VectorXd x_eval(nvars);
	for (Eigen::Index d = 0; d < nvars; ++d) {
		x_eval[d] = x[d];
	}

	double fval = NA_REAL;
	Eigen::VectorXd grad(nvars);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();
	{
		cppad_parallel_setup(static_cast<std::size_t>(num_threads));
		funObj.get_fdfh(x_eval, fval, grad, H);
	}

	InnerOptResult out;
	out.fval = fval;
	out.gradient = grad;
	out.hessian = H;
	return out;
}

//' @rdname innerOpt
//' @export
// [[Rcpp::export]]
Rcpp::List all_derivs(
	const Rcpp::NumericVector& x,
	const Rcpp::List& ad_fun,
	const Rcpp::List& config)
{
	const Config config_c(config);
	ad_groups* groups_ptr = get_ad_groups(ad_fun);
	std::vector<double> x_vec(x.begin(), x.end());

	const InnerOptResult result_c = all_derivs(x_vec, *groups_ptr, config_c);

	return Rcpp::List::create(
		Rcpp::Named("fval") = Rcpp::wrap(result_c.fval),
		Rcpp::Named("gradient") = Rcpp::wrap(result_c.gradient),
		Rcpp::Named("hessian") = eigen_to_dgCMatrix(result_c.hessian)
	);
}

InnerOptResult inner_opt(
	const std::vector<double>& parameters,
	const std::vector<double>& gamma,
	ad_groups& groups,
	const Config& config,
	const TrustControl& control,
	bool deriv)
{
	using Tvec = Eigen::VectorXd;
	using THess = Eigen::SparseMatrix<double>;
	using TPreLLt = Eigen::SimplicialLLT<THess>;

	if (parameters.size() != config.Nbeta + config.Ntheta) {
		Rcpp::stop(
			"parameters has length %d but expected Nbeta+Ntheta=%d",
			static_cast<int>(parameters.size()),
			static_cast<int>(config.Nbeta + config.Ntheta)
		);
	}
	if (gamma.size() != config.Ngamma) {
		Rcpp::stop(
			"gamma has length %d but expected Ngamma=%d",
			static_cast<int>(gamma.size()),
			static_cast<int>(config.Ngamma)
		);
	}
	if (config.Ngamma == 0) {
		Rcpp::stop("Ngamma must be > 0");
	}

	const int num_threads = config.num_threads > 0 ? config.num_threads : 1;

	Tvec gamma_start(config.Ngamma);
	Tvec solution(config.Ngamma);
	Tvec fullParams(config.Nparams);
	Tvec grad(config.Ngamma);

	std::copy(gamma.begin(), gamma.end(), gamma_start.data());

	std::vector<double> params_init(config.Nparams);
	for (size_t d = 0; d < config.Nbeta; ++d) {
		params_init[config.beta_begin + d] = parameters[d];
		fullParams[config.beta_begin + d] = parameters[d];
	}
	for (size_t d = 0; d < config.Ngamma; ++d) {
		params_init[config.gamma_begin + d] = gamma[d];
	}
	for (size_t d = 0; d < config.Ntheta; ++d) {
		params_init[config.theta_begin + d] = parameters[config.Nbeta + d];
		fullParams[config.theta_begin + d] = parameters[config.Nbeta + d];
	}

	AD_Func_Opt funObj(groups, params_init, true, num_threads);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();

	AD_Func_Opt funObjOuter(groups, params_init, false, num_threads);
	Eigen::SparseMatrix<double> Houter;
	Tvec gradOuter;

	double fval = NA_REAL;
	double radius = NA_REAL;
	int iterations = NA_INTEGER;
	MB_Status status = SUCCESS;

	InnerOptResult out;
	out.x_out.assign(groups.chol_pattern.i.size(), 0.0);
	out.d_out.assign(H.rows(), 0.0);


	{
		cppad_parallel_setup(static_cast<std::size_t>(num_threads));
		
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
		if (status != SUCCESS || grad_l2_sq > 10.0) {
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

		for (size_t d = 0; d < config.Ngamma; ++d) {
			fullParams[config.gamma_begin + d] = solution[d];
		}

		if (deriv) {
			gradOuter = Tvec(config.Nparams);
			Houter = funObjOuter.Htemplate.cast<double>();
			funObjOuter.get_fdfh(fullParams, fval, gradOuter, Houter);
		}
	}


	H.makeCompressed();
	const double log_det = chol_update(
		H,
		groups.chol_pattern.perm,
		groups.chol_pattern.p,
		groups.chol_pattern.i,
		out.x_out,
		out.d_out
	);

	out.fval = fval;
	out.solution = solution;
	out.full_parameters = fullParams;
	out.iterations = iterations;
	out.status = status;
	out.trust_radius = radius;

	if (R_finite(log_det)) {
		out.half_log_det = 0.5 * log_det;
		out.log_lik = -out.fval + config.Ngamma * ONEHALFLOGTWOPI - out.half_log_det;
		out.chol_ok = true;
	} else {
		Rcpp::warning(
			"inner_opt: non-positive or non-finite D; half_log_det not computed"
		);
	}

	if (deriv) {
		out.gradient = gradOuter;
		out.hessian = Houter;
	} else {
		out.gradient = grad;
		out.hessian = H;
	}

	return out;
}


static bool resolve_deriv(SEXP deriv_arg, const Rcpp::List& config)
{
	if (!Rf_isNull(deriv_arg)) {
		return Rcpp::as<bool>(deriv_arg);
	}
	if (config.containsElementNamed("deriv") && !Rf_isNull(config["deriv"])) {
		return Rcpp::as<bool>(config["deriv"]);
	}
	if (config.containsElementNamed("inner_only") && !Rf_isNull(config["inner_only"])) {
		return !Rcpp::as<bool>(config["inner_only"]);
	}
	return true;
}

//' @rdname innerOpt
//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
	const Rcpp::NumericVector parameters,
	const Rcpp::NumericVector gamma,
	const Rcpp::List& config,
	const Rcpp::List& ad_fun,
	SEXP control = R_NilValue,
	SEXP deriv = R_NilValue)
{
	try {
		if (ad_fun.size() == 0) {
			Rcpp::stop("inner_opt requires a non-empty ad_fun list");
		}

		const Config config_c(config);
		const Rcpp::List control_list = (control == R_NilValue)
			? Rcpp::List()
			: Rcpp::as<Rcpp::List>(control);
		const TrustControl control_c(control_list);
		ad_groups* groups_ptr = get_ad_groups(ad_fun);

		const bool deriv_flag = resolve_deriv(deriv, config);

		std::vector<double> parameters_vec(parameters.begin(), parameters.end());
		std::vector<double> gamma_vec(gamma.begin(), gamma.end());

		const InnerOptResult result = inner_opt(
			parameters_vec,
			gamma_vec,
			*groups_ptr,
			config_c,
			control_c,
			deriv_flag
		);

		const CholPattern& chol_pat = groups_ptr->chol_pattern;

		return Rcpp::List::create(
			Rcpp::Named("log_lik") = Rcpp::wrap(result.log_lik),
			Rcpp::Named("fval") = Rcpp::wrap(result.fval),
			Rcpp::Named("solution") = Rcpp::wrap(result.solution),
			Rcpp::Named("full_parameters") = Rcpp::wrap(result.full_parameters),
			Rcpp::Named("gradient") = Rcpp::wrap(result.gradient),
			Rcpp::Named("hessian") = eigen_to_dgCMatrix(result.hessian),
			Rcpp::Named("chol_inner") = chol_inner_list(
				chol_pat,
				result.x_out,
				result.d_out
			),
			Rcpp::Named("half_log_det") = Rcpp::wrap(result.half_log_det),
			Rcpp::Named("iterations") = Rcpp::wrap(result.iterations),
			Rcpp::Named("status") = Rcpp::wrap(std::string(MB_strerror(result.status))),
			Rcpp::Named("trust.radius") = Rcpp::wrap(result.trust_radius),
			Rcpp::Named("method") = Rcpp::wrap("Sparse")
		);
	}
	catch (const Rcpp::exception& e) {
		Rcpp::stop("inner_opt failed: %s", e.what());
	}
	catch (...) {
		Rcpp::stop("inner_opt failed with unknown error");
	}
}
