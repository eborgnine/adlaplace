// Adapted from trustOptim's CG-sparse.h (Michael Braun, MPL 2.0).
// The stock modified-Cholesky preconditioner equilibrates the Hessian and
// then factors that scaled matrix, while conjugate gradients call
// solve() as if the factor belonged to the unscaled Hessian. The search
// direction is then D H^{-1} D g rather than H^{-1} g, and CG stalls when
// column norms differ. This copy factors the Hessian itself.

#ifndef ADLAPLACE_CG_SPARSE_TRUSTOPTIM_H
#define ADLAPLACE_CG_SPARSE_TRUSTOPTIM_H

#include <CG-base.h>

using Eigen::Dynamic;

template<typename TP, typename TFunc, typename THess, typename TPreLLt>
class Trust_CG_Sparse : public Trust_CG_Base<TP, TFunc, THess, TPreLLt> {

    typedef typename TP::Index Index;
    typedef Trust_CG_Base<TP, TFunc, THess, TPreLLt> BaseType;

    typedef SparseMatrix<double> SparseMatrixXd;

    using BaseType::Bk;
    using BaseType::nvars;
    using BaseType::func;
    using BaseType::function_scale_factor;
    using BaseType::xk;
    using BaseType::precond_ID;
    using BaseType::PrecondLLt;
    using BaseType::status;
    using BaseType::rad;
    using BaseType::iter;

private:

    Index nnz;

    void update_precond() {

	switch (precond_ID) {
	case 0:
	    update_precond_identity();
	    break;
	case 1:
	    update_precond_modified_Cholesky();
	    break;
	default:
	    update_precond_identity();
	}
	return;
    }

    void update_precond_identity() {
	return;
    }

    void update_precond_modified_Cholesky() {

	bool success = false;
	double alpha = 0, beta, bmin;

	SparseMatrixXd BB(nvars, nvars);
	BB = Bk.template selfadjointView<Lower>();

	beta = sqrt(BB.cwiseAbs2().sum());
	bmin = BB.coeff(0, 0);
	for (int j = 0; j < nvars; j++) {
	    bmin = std::min(bmin, BB.coeff(j, j));
	}

	if (!(bmin > 0)) alpha = beta / 2.;

	for (int ii = 0; ii < 40 && !success; ++ii) {
	    PrecondLLt.factorize(BB);
	    if (PrecondLLt.info() == Eigen::Success) {
		success = true;
	    } else {
		alpha = std::max(2 * alpha, beta / 2.) - alpha;
		for (int j = 0; j < nvars; j++) {
		    BB.coeffRef(j, j) += alpha;
		}
	    }
	}
    }

    void init_precond() {
	switch (precond_ID) {
	case 0:
	    init_precond_identity();
	    break;
	case 1:
	    init_precond_modified_Cholesky();
	    break;
	default:
	    init_precond_identity();
	}
	return;
    }

    void init_precond_identity() {
	SparseMatrixXd tmp(nvars, nvars);
	tmp.reserve(nvars);

	for (int j = 0; j < nvars; j++) {
	    tmp.startVec(j);
	    tmp.insertBack(j, j) = 1.;
	}
	tmp.makeCompressed();

	PrecondLLt.compute(tmp);

	return;
    }

    inline void init_precond_modified_Cholesky() {
	PrecondLLt.analyzePattern(Bk);
	update_precond_modified_Cholesky();
    }

    void init_sparse_hessian() {
	Bk.resize(nvars, nvars);
	Bk.reserve(nnz);
	func.get_hessian(xk, Bk);
    }

    inline void update_sparse_hessian() {
	func.get_hessian(xk, Bk);
    }

    inline void update_hessian() {
	update_sparse_hessian();
	Bk *= function_scale_factor;
    }

public:
    Trust_CG_Sparse(TFunc & func_,
		    const MatrixBase<TP>& startX_,
		    const double & rad_,
		    const double & min_rad_,
		    const double & tol_,
		    const double & prec_,
		    const int & report_freq_,
		    const int  report_level_,
		    const int & header_freq_,
		    const int & report_precision_,
		    const int & maxit_,
		    const double & contract_factor_,
		    const double & expand_factor_,
		    const double & contract_threshold_,
		    const double & expand_threshold_rad_,
		    const double & expand_threshold_ap_,
		    const double & function_scale_factor_,
		    const int & precond_refresh_freq_,
		    const int & precond_ID_,
		    const int & trust_iter_) :
	BaseType(func_, startX_, rad_, min_rad_, tol_, prec_, report_freq_,
		 report_level_, header_freq_, report_precision_, maxit_, contract_factor_,
		 expand_factor_, contract_threshold_, expand_threshold_rad_,
		 expand_threshold_ap_, function_scale_factor_,
		 precond_refresh_freq_, precond_ID_, trust_iter_),
	nnz(func_.get_nnz())
    {
	init_sparse_hessian();
	Bk *= function_scale_factor;
	init_precond();
    }

    template<typename Tvec, typename Tout>
    MB_Status get_current_state(const MatrixBase<Tvec> & pars_,
				double & fval,
				const MatrixBase<Tvec> & grad_,
				const SparseMatrixBase<Tout> & hess_,
				int & iterations,
				double & radius) {

	MatrixBase<Tvec> & pars = const_cast<MatrixBase<Tvec>& >(pars_);
	MatrixBase<Tvec> & grad = const_cast<MatrixBase<Tvec>& >(grad_);
	SparseMatrixBase<Tout> & hess = const_cast<SparseMatrixBase<Tout>&>(hess_);

	pars = xk;
	func.get_fdf(xk, fval, grad);
	func.get_hessian(xk, hess);

	iterations = iter;
	radius = rad;

	return status;
    }

}; // end class Trust_CG_Sparse

#endif
