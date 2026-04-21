
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
//'   \code{gradient} (full/outer gradient at the optimized \code{gamma}),
//'   \code{hessian} (sparse list), \code{iterations}, \code{status},
//'   \code{trust.radius}, and \code{method}.}
//'
//' @details
//' This calls the sparse method from the \code{TrustOptim} package via the Cpp interface.  
//'
//' @name innerOpt


#include <Rcpp.h>
#include <Eigen/Sparse>
#include <cmath>
#include <cstdio>
#include <cstdlib>

// from trustOptim
#include <CG-sparse.h> 

#include "adlaplace/runtime/rviews.hpp"
#include "adlaplace/ompad.hpp"
#include "adlaplace/math/constants.hpp"
#include "trustOptimWrappers.hpp"

auto get_double_ctrl = [](const Rcpp::List& ctl, const char* key, double def) {
	return ctl.containsElementNamed(key) ? Rcpp::as<double>(ctl[key]) : def;
};
auto get_int_ctrl = [](const Rcpp::List& ctl, const char* key, int def) {
	return ctl.containsElementNamed(key) ? Rcpp::as<int>(ctl[key]) : def;
};

Rcpp::List eigen_to_list(
	const Eigen::SparseMatrix<double> &M,
	const bool upper = true) {

	Eigen::SparseMatrix<double> Mcopy;
	if (upper) {
		Mcopy = M.triangularView<Eigen::Upper>();
	} else {
		Mcopy = M;
	}
	Mcopy.makeCompressed();


	const Eigen::Index nnz  = Mcopy.nonZeros();
	const Eigen::Index ncolP1 = Mcopy.cols()+1;


	Rcpp::NumericVector x(nnz);
	Rcpp::IntegerVector i(nnz);
	Rcpp::IntegerVector p(ncolP1);

    // copy from Eigen into R vectors
	std::copy(Mcopy.valuePtr(),
		Mcopy.valuePtr() + nnz,
		x.begin());

	std::copy(Mcopy.innerIndexPtr(),
		Mcopy.innerIndexPtr() + nnz,
		i.begin());

	std::copy(Mcopy.outerIndexPtr(),
		Mcopy.outerIndexPtr() + ncolP1,
		p.begin());

	Rcpp::List result = Rcpp::List::create(
		Rcpp::_["i"] = i, 
		Rcpp::_["p"] = p,
		Rcpp::_["x"] = x,
		Rcpp::_["dims"] = Rcpp::IntegerVector::create(
			static_cast<int>(M.rows()),
			static_cast<int>(M.cols())
			),
		Rcpp::_["index1"] = false
		);
	if (upper) {
		result["symmetric"] = true;
	} else {
		result["triangular"] = true;
	}
	return(result);
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
		Rcpp::Named("fval")          = Rcpp::wrap(fval),
		Rcpp::Named("gradient")      = Rcpp::wrap(grad),
		Rcpp::Named("hessian")       = eigen_to_list(H)
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


	using Tvec   = Eigen::VectorXd;
	using THess   = Eigen::SparseMatrix<double>; 
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

	Tvec gamma_start(configC.Ngamma), solution(configC.Ngamma);
	Tvec fullParams(configC.Nparams);
	Tvec grad(configC.Ngamma), gradOuter(configC.Nparams);

	for (size_t d = 0; d < configC.Ngamma; ++d) {
		gamma_start[d] = gamma[d];
	}

	std::vector<double> params_init(configC.Nparams);
	// copy in beta and theta from parameters, and gamma from the gamma argument
	for (size_t d = 0; d < configC.Nbeta; ++d) {
		params_init[configC.beta_begin + d] = parameters[d];
		fullParams[configC.beta_begin + d] = parameters[d];
	}
	for (size_t d = 0; d < configC.Ngamma; ++d) {
		params_init[configC.gamma_begin + d] = gamma[d];
	}
	for (size_t d = 0; d < configC.Ntheta; ++d) {
		params_init[configC.theta_begin + d] = parameters[configC.Nbeta + d];
		fullParams[configC.theta_begin + d] = parameters[configC.Nbeta +d];
	}



	AD_Func_Opt 
	funObj(
		adFun,
		params_init,
		true,
		num_threads),
	funObjOuter(
		adFun,
		params_init,
		false,           // inner=false
		num_threads);


	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>(),
		Houter = funObjOuter.Htemplate.cast<double>();


	double fval = NA_REAL, radius = NA_REAL;
	int iterations = NA_INTEGER;
	MB_Status status;

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
		// get full hessian, gradient
		funObjOuter.get_fdfh(fullParams, fval, gradOuter, Houter);
	}

	Rcpp::NumericVector solutionR(solution.size());
	for(size_t D=0;D<solution.size();D++) {
		solutionR[D] = solution[D];
	}
	Rcpp::NumericVector gradientR(gradOuter.size());
	for(size_t D=0;D<gradientR.size();D++) {
		gradientR[D] = gradOuter[D];
	}

	Rcpp::List res = Rcpp::List::create(
		Rcpp::Named("fval") = Rcpp::wrap(fval),
		Rcpp::Named("solution") = solutionR,
		Rcpp::Named("gradient") = gradientR,
		Rcpp::Named("hessian") = eigen_to_list(Houter),
		Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		Rcpp::Named("method") = Rcpp::wrap("Sparse")
		);

	return(res);
}
