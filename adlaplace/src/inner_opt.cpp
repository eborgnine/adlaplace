
//' Inner optimization over gamma using trust-region CG (sparse)
//'
//' Runs the inner optimization problem (typically over \eqn{\gamma}) using the
//' trustOptim sparse trust-region Conjugate Gradient solver. This function
//' evaluates the objective, gradient, and Hessian through the pre-built AD pack
//' (external pointer) and returns the solution along with curvature information.
//'
//' @param x Numeric full parameter vector of length \code{Nparams}, used by
//'   \code{all_derivs()} for outer/full derivatives.
//' @param parameters Numeric vector of fixed outer parameters
//'   (\code{beta}, \code{theta}; length \code{Nbeta+Ntheta}) used by
//'   \code{inner_opt()}.
//' @param gamma Numeric vector of starting values for inner parameters
//'   (\code{gamma}; length \code{Ngamma}) used by \code{inner_opt()}.
//' @param data Data list used to construct the AD backend when \code{adPack}
//'   is not supplied.
//' @param adPack External pointer created by \code{getAdFun()} (class
//'   \code{"adlaplace_handle_ptr"}), or a list containing element \code{adFun}.
//' @param config Configuration list. Must include \code{gamma}, fixed
//'   \code{beta}/\code{theta}, and group/sparsity settings.
//' @param control List of trust-region control parameters (see
//'   \code{trustOptim}).
//' @param adPack Optional backend handle/list from \code{getAdFun()}.
//'   If provided, it is reused.
//'
//' @return A list with components:
//' \itemize{
//'   \item \code{all_derivs}: list with \code{fval}, \code{gradient}, and
//'         \code{hessian} for the outer/full derivatives (\code{inner=FALSE}).
//'   \item \code{minusLogLik}: scalar \eqn{-\ell(\hat\gamma)} plus the Laplace
//'         correction \eqn{\tfrac{1}{2}\log|H| + \tfrac{n}{2}\log(2\pi)}.
//'   \item \code{fval}: scalar objective at the solution (typically \eqn{-\ell}).
//'   \item \code{halfLogDet}: \eqn{\tfrac{1}{2}\log|H|} from sparse LDLT.
//'   \item \code{solution}: optimized parameter vector (length \code{Ngamma}).
//'   \item \code{gradient}: gradient at solution (length \code{Ngamma}).
//'   \item \code{hessian}: Hessian as a dgCMatrix-like list with slots
//'         \code{i,p,x,Dim} (0-based indices).
//'   \item \code{cholHessian}: sparse LDLT factors as a list with
//'         \code{P} (permutation indices), \code{D} (diagonal), and
//'         \code{L} (lower-triangular factor in dgCMatrix-like form).
//'   \item \code{iterations}: number of trust-region iterations.
//'   \item \code{status}: solver status string.
//'   \item \code{trust.radius}: final trust-region radius.
//'   \item \code{method}: character, here \code{"Sparse"}.
//' }
//'
//' @details
//' This calls the sparse method from the \code{TrustOptim} package via the Cpp interface.  
//'
//' @name innerOpt


#include <Rcpp.h>
#include <Eigen/Sparse>

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
	SEXP adPack,
	const Rcpp::List& config)
{

	const Config configC(config);
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

	AD_Func_Opt funObj(adPack, params_init, false); // inner=false
	const Eigen::Index nvars = static_cast<Eigen::Index>(funObj.get_nvars());
	Eigen::VectorXd x_eval(nvars);
	for (Eigen::Index d = 0; d < nvars; ++d) {
		x_eval[d] = x[d];
	}

	double fval = NA_REAL;
	Eigen::VectorXd grad(nvars);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();
	
	{
		auto guard=cppad_parallel_setup(configC.num_threads);
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
	SEXP adPack = R_NilValue
	) {
	const Config configC(config);
	if (adPack == R_NilValue) {
		Rcpp::stop("inner_opt requires a non-NULL adPack");
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

	Tvec gamma_start(configC.Ngamma);
	for (size_t d = 0; d < configC.Ngamma; ++d) {
		gamma_start[d] = gamma[d];
	}

	std::vector<double> params_init(configC.Nparams);
	// copy in beta and theta from parameters, and gamma from the gamma argument
	for (size_t d = 0; d < configC.Nbeta; ++d) {
		params_init[configC.beta_begin + d] = parameters[d];
	}
	for (size_t d = 0; d < configC.Ngamma; ++d) {
		params_init[configC.gamma_begin + d] = gamma[d];
	}
	for (size_t d = 0; d < configC.Ntheta; ++d) {
		params_init[configC.theta_begin + d] = parameters[configC.Nbeta + d];
	}

	AD_Func_Opt funObj(adPack, params_init, true);


	Trust_CG_Sparse<Tvec, AD_Func_Opt, THess, TPreLLt> opt(
		funObj, 
		gamma_start,
		rad,
		min_rad,
		tol,
		prec,
		report_freq,
		report_level,
		header_freq,
		report_precision,
		maxit,
		contract_factor,
		expand_factor,
		contract_threshold,
		expand_threshold_rad,
		expand_threshold_ap,
		function_scale_factor,
		precond_refresh_freq,
		precond_ID,
//		ctrl.quasi_newton_method,
		trust_iter
		);


	if(configC.verbose) {
		Rcpp::Rcout << "starting opt..";
	}

	{
		auto guard=cppad_parallel_setup(configC.num_threads);
		opt.run();
	}

	if(configC.verbose) {
		Rcpp::Rcout << ".done.\n";
	}


	Tvec P(configC.Ngamma);
	Tvec grad(configC.Ngamma);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();

	double fval = NA_REAL, radius = NA_REAL;
	int iterations = NA_INTEGER;
	MB_Status status;


	status = opt.get_current_state(P, fval, grad, H, iterations, radius);

	// get log determinant of hessian
	Eigen::SparseMatrix<double> Htri = H.triangularView<Eigen::Upper>();
	Htri.makeCompressed();
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper> ldlt;
	ldlt.compute(Htri);
	if (ldlt.info() != Eigen::Success) {
		Rcpp::warning("LDLT factorization failed; H may not be SPD");
	}

	// log likelhood

	const auto& D = ldlt.vectorD();
	double logdetH = 0.0;
	for (int k = 0; k < D.size(); ++k) {
		logdetH += std::log(D[k]);
	}
	const double halfLogDet = logdetH/2;

	const double minusLogLik = fval + halfLogDet + configC.Ngamma * ONEHALFLOGTWOPI;  

 	 // save chol
	const auto& Pidx = ldlt.permutationP().indices();    
	const Eigen::SparseMatrix<double> L = ldlt.matrixL();

	Rcpp::NumericVector D_R(D.size());
	for (int k = 0; k < D.size(); ++k) {
		D_R[k] = D[k];
	}
	Rcpp::IntegerVector P_R(Pidx.size());
	for (int k = 0; k < Pidx.size(); ++k){
		P_R[k] = Pidx[k];
	}

	Rcpp::List cholHessian = Rcpp::List::create(
		Rcpp::Named("P") = P_R,
		Rcpp::Named("D") = D_R,
		Rcpp::Named("L") = eigen_to_list(L, false)
		);

	Rcpp::NumericVector solution(P.size());
	for(size_t D=0;D<solution.size();D++) {
		solution[D] = P[D];
	}

	Rcpp::List res = Rcpp::List::create(
		Rcpp::Named("minusLogLik") = Rcpp::wrap(minusLogLik),
		Rcpp::Named("fval") = Rcpp::wrap(fval),
		Rcpp::Named("halfLogDet") = Rcpp::wrap(halfLogDet),
		Rcpp::Named("solution") = Rcpp::wrap(solution),
		Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		Rcpp::Named("hessian") = eigen_to_list(Htri),
		Rcpp::Named("cholHessian") = cholHessian,
		Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		Rcpp::Named("method") = Rcpp::wrap("Sparse")
		);

	return(res);

}
