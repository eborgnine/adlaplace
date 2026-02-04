
//' Inner optimization over gamma using trust-region CG (sparse)
//'
//' Runs the inner optimization problem (typically over \eqn{\gamma}) using the
//' trustOptim sparse trust-region Conjugate Gradient solver. This function
//' evaluates the objective, gradient, and Hessian through the pre-built AD pack
//' (external pointer) and returns the solution along with curvature information.
//'
//' @param start Numeric vector of starting values for the inner parameters
//'   (length \code{Ngamma}).
//' @param adPack External pointer created by \code{getAdFun()} (class
//'   \code{"adpack_ptr"}). It contains per-group AD tapes and workspaces.
//' @param config List of configuratio.  Must include  \code{gamma} (starting values), fixed values of \code{beta} and
//'   \code{theta}, and sparsity/match info under \code{config$sparsity}.
//' @param control List of trust-region control parameters (see
//'   \code{trustOptim}).
//'
//' @return A list with components:
//' \itemize{
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

#include "adlaplace/utils.hpp"
#include "adlaplace/foromp.hpp"
#include "adlaplace/constants.hpp"
#include "trustOptimWrappers.hpp"
#include "trustOptimControl.hpp"

Rcpp::List eigen_to_list(
	const Eigen::SparseMatrix<double> &M) {

	const Eigen::Index nnz  = M.nonZeros();
	const Eigen::Index ncolP1 = M.cols()+1;

	Rcpp::NumericVector x(nnz);
	Rcpp::IntegerVector i(nnz);
	Rcpp::IntegerVector p(ncolP1);

    // copy from Eigen into R vectors
	std::copy(M.valuePtr(),
		M.valuePtr() + nnz,
		x.begin());

	std::copy(M.innerIndexPtr(),
		M.innerIndexPtr() + nnz,
		i.begin());

	std::copy(M.outerIndexPtr(),
		M.outerIndexPtr() + ncolP1,
		p.begin());

	Rcpp::List result = Rcpp::List::create(
		Rcpp::_["i"] = i, 
		Rcpp::_["p"] = p,
		Rcpp::_["x"] = x
		);
	return(result);
}



//' @rdname innerOpt
//' @export
// [[Rcpp::export]]
Rcpp::List innerOptTest(
	SEXP adPack,
	const Rcpp::List& config
	)
{

	const Config configC(config);
	const size_t Nparams = configC.Ngamma;

	Eigen::VectorXd parameters(Nparams); 
	for(size_t D=0;D<Nparams;D++) {
		parameters[D] = configC.gamma[D];
	}

	AD_Func_Opt funObj(adPack, configC);

	double fval = NA_REAL;
	Eigen::VectorXd grad(Nparams);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();
	
	{
		auto guard=cppad_parallel_setup(configC.num_threads);
		funObj.get_fdfh(parameters, fval, grad, H);
	}


  // mirror innerOpt list structure as closely as possible
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
Rcpp::List innerOpt(
	SEXP adPack,
	const Rcpp::List& config,
	const Rcpp::List& control
	)
{

	const TrustControl ctrl(control); 
	const Config configC(config);

	using Tvec   = Eigen::VectorXd;
	using THess   = Eigen::SparseMatrix<double>; 
	using TPreLLt = Eigen::SimplicialLLT<THess>;

	const size_t Nparams = configC.Ngamma;


	Tvec parameters(Nparams); 
	for(size_t D=0;D<Nparams;D++) {
		parameters[D] = configC.gamma[D];
	}

	AD_Func_Opt funObj(adPack, configC);


	Trust_CG_Sparse<Tvec, AD_Func_Opt, THess, TPreLLt> opt(
		funObj, 
		parameters,
		ctrl.rad,
		ctrl.min_rad,
		ctrl.tol,
		ctrl.prec,
		ctrl.report_freq,
		ctrl.report_level,
		ctrl.header_freq,
		ctrl.report_precision,
		ctrl.maxit,
		ctrl.contract_factor,
		ctrl.expand_factor,
		ctrl.contract_threshold,
		ctrl.expand_threshold_rad,
		ctrl.expand_threshold_ap,
		ctrl.function_scale_factor,
		ctrl.precond_refresh_freq,
		ctrl.precond_ID,
//		ctrl.quasi_newton_method,
		ctrl.trust_iter
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


	Tvec P(Nparams);
	Tvec grad(Nparams);
	Eigen::SparseMatrix<double> H = funObj.Htemplate.cast<double>();

	double fval = NA_REAL, radius = NA_REAL;
	int iterations = NA_INTEGER;
	MB_Status status;


	status = opt.get_current_state(P, fval, grad, H, iterations, radius);

	// get log determinant of hessian
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
	ldlt.compute(H);
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

	const double minusLogLik = fval + halfLogDet + Nparams * ONEHALFLOGTWOPI;  

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
		Rcpp::Named("L") = eigen_to_list(L)
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
		Rcpp::Named("hessian") = eigen_to_list(H),
		Rcpp::Named("cholHessian") = cholHessian,
		Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		Rcpp::Named("method") = Rcpp::wrap("Sparse")
		);

	return(res);

}

