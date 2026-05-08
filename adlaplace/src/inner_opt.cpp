
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
//' @param adFun External pointer returned by \code{getAdFun()}.
//' @param config Configuration list with model dimensions, groups, and
//'   sparsity information.
//' @param control List of trust-region control parameters for
//'   \code{inner_opt()} (see \pkg{trustOptim}).
//'
//' @return
//'   \item{\code{all_derivs()}}{Returns \code{fval}, \code{gradient}, and
//'   \code{hessian} for full/outer derivatives at \code{x}.}
//'   \item{\code{inner_opt()}}{Returns \code{fval}, \code{solution},
//'   \code{gradient} (inner gradient if \code{inner_only} is set, otherwise full gradient),
//'   \code{hessian} (inner Hessian if \code{inner_only} is set, otherwise full Hessian),
//'   \code{iterations}, \code{status}, \code{trust.radius}, and \code{method}.}
//'
//' @details
//' This calls the sparse method from the \code{TrustOptim} package via the Cpp interface.  
//'
//' @name innerOpt

// Standard C
#include <cmath>
#include <cstdio>
#include <cstdlib>

// Eigen
#include <Eigen/Sparse>

// Rcpp
#include <Rcpp.h>

// trustOptim
#include <CG-sparse.h>

// Local
#include "adlaplace/math/constants.hpp"
#include "adlaplace/runtime/rviews.hpp"
#include "adlaplace/ompad.hpp"
#include "trustOptimWrappers.hpp"

auto get_double_ctrl = [](const Rcpp::List& ctl, const char* key, double def) {
	if (ctl.containsElementNamed(key) && !Rf_isNull(ctl[key])) {
		return Rcpp::as<double>(ctl[key]);
	}
	return def;
	};
auto get_int_ctrl = [](const Rcpp::List& ctl, const char* key, int def) {
	if (ctl.containsElementNamed(key) && !Rf_isNull(ctl[key])) {
		return Rcpp::as<int>(ctl[key]);
	}
	return def;
	};



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


//' @rdname innerOpt
//' @export
// [[Rcpp::export]]
Rcpp::List all_derivs(
	const Rcpp::NumericVector& x,
	SEXP adFun,
	const Rcpp::List& config)
{

	const Config configC(config);
	const int num_threads = configC.num_threads > 0 ? configC.num_threads : 1;
	if (x.size() != static_cast<R_xlen_t>(configC.Nparams)) {
		Rcpp::stop(
			"all_derivs: x has length %d but expected Nparams=%d",
			static_cast<int>(x.size()),
			static_cast<int>(configC.Nparams)
		);
	}
	std::vector<double> params_init(configC.Nparams);
	for (size_t d = 0; d < configC.Nparams; ++d) {
		params_init[d] = x[d];
	}

	AD_Func_Opt funObj(
		adFun,
		params_init,
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


	// mirror inner_opt list structure as closely as possible
	Rcpp::List res = Rcpp::List::create(
		Rcpp::Named("fval") = Rcpp::wrap(fval),
		Rcpp::Named("gradient") = Rcpp::wrap(grad),
		Rcpp::Named("hessian") = eigen_to_dgCMatrix(H)
	);

	return res;
}



//' @rdname innerOpt
//' @export
// [[Rcpp::export]]
Rcpp::List inner_opt(
	const Rcpp::NumericVector parameters, // beta and theta, fixed
	const Rcpp::NumericVector gamma, // starting values
	const Rcpp::List& config,
	const Rcpp::List& control,
	SEXP adFun = R_NilValue
) {
	try {
		const Config configC(config);
	if (adFun == R_NilValue) {
		Rcpp::stop("inner_opt requires a non-NULL adFun");
	}

	const double rad = get_double_ctrl(control, "step.size", 1.0);
	const double min_rad = get_double_ctrl(control, "min.step.size", 1e-8);
	const double tol = get_double_ctrl(control, "cg.tol", 1e-4);
	const double prec = get_double_ctrl(control, "grad.tol", 1e-6);
	const int report_freq = get_int_ctrl(control, "report.freq", 0);
	const int report_level = get_int_ctrl(control, "report.level", 0);
	const int header_freq = get_int_ctrl(control, "header.freq", 10);
	const int report_precision = get_int_ctrl(control, "report.precision", 6);
	const int maxit = get_int_ctrl(control, "maxit", 100);
	const double contract_factor = get_double_ctrl(control, "contract.factor", 0.5);
	const double expand_factor = get_double_ctrl(control, "expand.factor", 2.0);
	const double contract_threshold = get_double_ctrl(control, "contract.threshold", 0.25);
	const double expand_threshold_rad = get_double_ctrl(control, "expand.threshold.rad", 0.8);
	const double expand_threshold_ap = get_double_ctrl(control, "expand.threshold.ap", 0.75);
	const double function_scale_factor = get_double_ctrl(control, "function.scale.factor", 1.0);
	const int precond_refresh_freq = get_int_ctrl(control, "precond.refresh", 5);
	const int precond_ID = get_int_ctrl(control, "precond.ID", 0);
	const int trust_iter = get_int_ctrl(control, "trust.iter", 50);
	const int num_threads = configC.num_threads > 0 ? configC.num_threads : 1;


	using Tvec = Eigen::VectorXd;
	using THess = Eigen::SparseMatrix<double>;
	using TPreLLt = Eigen::SimplicialLLT<THess>;

	if (parameters.size() != static_cast<R_xlen_t>(configC.Nbeta + configC.Ntheta)) {
		Rcpp::stop(
			"parameters has length %d but expected Nbeta+Ntheta=%d",
			static_cast<int>(parameters.size()),
			static_cast<int>(configC.Nbeta + configC.Ntheta)
		);
	}
	if (gamma.size() != static_cast<R_xlen_t>(configC.Ngamma)) {
		Rcpp::stop(
			"gamma has length %d but expected Ngamma=%d",
			static_cast<int>(gamma.size()),
			static_cast<int>(configC.Ngamma)
		);
	}
	if (configC.Ngamma == 0) {
		Rcpp::stop("Ngamma must be > 0");
	}

	Tvec gamma_start(configC.Ngamma), solution(configC.Ngamma);
	Tvec fullParams(configC.Nparams);
	Tvec grad(configC.Ngamma);

	// Copy gamma to gamma_start
	std::copy(gamma.begin(), gamma.end(), gamma_start.data());

	std::vector<double> params_init(configC.Nparams);
	// Copy beta, gamma, theta into params_init and fullParams
	for (size_t d = 0; d < configC.Nbeta; ++d) {
		params_init[configC.beta_begin + d] = parameters[d];
		fullParams[configC.beta_begin + d] = parameters[d];
	}
	for (size_t d = 0; d < configC.Ngamma; ++d) {
		params_init[configC.gamma_begin + d] = gamma[d];
	}
	for (size_t d = 0; d < configC.Ntheta; ++d) {
		params_init[configC.theta_begin + d] = parameters[configC.Nbeta + d];
		fullParams[configC.theta_begin + d] = parameters[configC.Nbeta + d];
	}

	AD_Func_Opt funObj(
		adFun,
		params_init,
		true,
		num_threads);

	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();

// cholesky template
    Rcpp::List adFun_list(adFun);
    Rcpp::List hessians = adFun_list["hessians"];
	bool use_chol_template;
	SEXP chol_template;
    if (hessians.containsElementNamed("chol_inner") && 
            !Rf_isNull(hessians["chol_inner"])) {
		use_chol_template = true;
		chol_template = hessians["chol_inner"];
	} else {
		use_chol_template = false;
	}

	AD_Func_Opt funObjOuter(
		adFun,
		params_init,
		false,           // inner=false
		num_threads);

	Eigen::SparseMatrix<double>	Houter;
	Tvec gradOuter;

	// Check if we should return inner-only results
	bool use_inner = false;
	if (config.containsElementNamed("inner_only") && !Rf_isNull(config["inner_only"])) {
		use_inner = Rcpp::as<bool>(config["inner_only"]);
	}

	double fval = NA_REAL, radius = NA_REAL;
	int iterations = NA_INTEGER;
	MB_Status status = SUCCESS;

	{
		cppad_parallel_setup(static_cast<std::size_t>(num_threads));
		Trust_CG_Sparse<Tvec, AD_Func_Opt, THess, TPreLLt> opt(
			funObj, gamma_start, rad, min_rad, tol, prec,
			report_freq, report_level, header_freq, report_precision,
			maxit, contract_factor, expand_factor, contract_threshold,
			expand_threshold_rad, expand_threshold_ap, function_scale_factor,
			precond_refresh_freq, precond_ID, trust_iter
		);

		opt.run();
		// H and grad won't be populated
		status = opt.get_current_state(solution, fval, grad,
			H, iterations, radius);

		// copy gamma to full parameters
		for (size_t d = 0; d < configC.Ngamma; ++d) {
			fullParams[configC.gamma_begin + d] = solution[d];
		}

		// get full hessian, gradient (only if not inner_only or inner_only is FALSE)
		if (!use_inner) {
			gradOuter = Tvec(configC.Nparams);
			Houter = funObjOuter.Htemplate.cast<double>();

			funObjOuter.get_fdfh(fullParams, fval, gradOuter, Houter);
		}

	} // parallel block

	// Compute Cholesky decomposition of inner Hessian using Matrix package

	static Rcpp::Function Matrix_Cholesky = Rcpp::Environment::namespace_env("Matrix")["Cholesky"];
	static Rcpp::Function Matrix_update = Rcpp::Environment::namespace_env("Matrix")["update"];
	static Rcpp::Function Matrix_determinant = Rcpp::Environment::namespace_env("Matrix")["determinant"];


	SEXP H_mat = eigen_to_dgCMatrix(H);
	SEXP H_chol;

	if(use_chol_template) {
		Rcpp::Rcout << "!";
		H_chol = Matrix_update(
			Rcpp::Named("object") = chol_template, 
			Rcpp::Named("parent") = H_mat);
	} else{
		H_chol = Matrix_Cholesky(H_mat, Rcpp::Named("perm") = true, Rcpp::Named("LDL") = true);
	}
    Rcpp::List log_det_list = Matrix_determinant(H_chol, Rcpp::Named("logarithm") = true, Rcpp::Named("sqrt") = true);

	double half_log_det = Rcpp::as<double>(log_det_list["modulus"]);
	double log_lik_val = -fval + configC.Ngamma * ONEHALFLOGTWOPI - half_log_det;



	// Conditionally return inner or outer gradient/Hessian
	// If config has "inner_only" and it's not FALSE, return inner H/grad
	SEXP hessian_out_mat;
	if (!use_inner) {
		hessian_out_mat = eigen_to_dgCMatrix(Houter);
	} else {
		hessian_out_mat = H_mat;
	}
	const Tvec& gradient_out = use_inner ? grad : gradOuter;

	Rcpp::NumericVector gradientR = Rcpp::wrap(gradient_out);
	Rcpp::NumericVector full_parameters_R = Rcpp::wrap(fullParams);
	Rcpp::NumericVector solutionR = Rcpp::wrap(solution);

	Rcpp::List res = Rcpp::List::create(
		Rcpp::Named("log_lik") = Rcpp::wrap(log_lik_val),
		Rcpp::Named("fval") = Rcpp::wrap(fval),
		Rcpp::Named("solution") = solutionR,
		Rcpp::Named("full_parameters") = full_parameters_R,
		Rcpp::Named("gradient") = gradientR,
		Rcpp::Named("hessian") = hessian_out_mat,
		Rcpp::Named("Hchol") = H_chol,
		Rcpp::Named("half_log_det") = Rcpp::wrap(half_log_det),
		Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		Rcpp::Named("status") = Rcpp::wrap(std::string(MB_strerror(status))),
		Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		Rcpp::Named("method") = Rcpp::wrap("Sparse")
	);

	return(res);
	} catch (const Rcpp::exception& e) {
		Rcpp::stop("inner_opt failed: %s", e.what());
	} catch (...) {
		Rcpp::stop("inner_opt failed with unknown error");
	}
}
