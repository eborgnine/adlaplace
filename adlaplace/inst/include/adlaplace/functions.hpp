#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "adlaplace/adfun.hpp"
#include <RcppEigen.h>

static const std::string JAC_COLOR = "cppad";  
static const std::string HESS_COLOR = "cppad.symmetric";


#ifdef DEBUG
#include<Rcpp.h>
#endif			


inline void grad(
	const CppAD::vector<double> &x,
	double &f,
	std::vector<double> &g, 
	std::vector<GroupPack>& adPack
	) {

	f = 0.0;
	const std::size_t n = static_cast<std::size_t>(x.size());
	g.assign(n, 0.0); 

  #pragma omp parallel 
	{
		const size_t Npack = adPack.size();
		double fOutHere = 0.0;
		std::vector<double> gradHere(n,0.0);
		CppAD::vector<double> w(1);
		w[0] = 1.0;

      #pragma omp for nowait
		for(int D=0;D<Npack;D++) {

			CppAD::vector<double> y = adPack[D].fun.Forward(0, x);
			fOutHere += y[0];
#ifdef DEBUG
			Rcpp::Rcout << "\n=== DEBUG: Group " << D << " ===\n";

    // Domain size of ADFun
			Rcpp::Rcout << "fun.Domain() = " << adPack[D].fun.Domain() << "\n";

    // Pattern info
			{
				const auto& pat = adPack[D].pattern_grad;
				Rcpp::Rcout << "pattern_grad: nr=" << pat.nr() 
				<< " nc=" << pat.nc() 
				<< " nnz=" << pat.nnz() << "\n";

				for (size_t k = 0; k < pat.nnz(); ++k) {
					Rcpp::Rcout << "  pattern[" << k << "]: row=" 
					<< pat.row()[k] 
					<< " col=" << pat.col()[k] << "\n";
				}
			}

    // out_grad info
			{
				const auto& out = adPack[D].out_grad;
				Rcpp::Rcout << "out_grad: nr=" << out.nr() 
				<< " nc=" << out.nc() 
				<< " nnz=" << out.nnz() << "\n";

				for (size_t k = 0; k < out.nnz(); ++k) {
					Rcpp::Rcout << "  out_grad[" << k 
					<< "]: row=" << out.row()[k]
					<< " col=" << out.col()[k]
					<< " val=" << out.val()[k] << "\n";
				}
			}

    // Check for index overflows
			{
				const auto& out = adPack[D].out_grad;
				for (size_t k = 0; k < out.nnz(); ++k) {
					if (out.col()[k] >= adPack[D].fun.Domain()) {
						Rcpp::Rcout << "!!! ERROR: out_grad.col[" << k 
						<< "] = " << out.col()[k] 
						<< " exceeds domain size!\n";
					}
				}
			}

			Rcpp::Rcout << "=== END DEBUG GROUP " << D << " ===\n";
#endif			


			adPack[D].fun.sparse_jac_rev(x,
				adPack[D].out_grad,   
				adPack[D].pattern_grad,
				JAC_COLOR,      
				adPack[D].work_grad);

#ifdef DEBUG
			Rcpp::Rcout << ".";
#endif			

			const auto& index = adPack[D].out_grad.col();    
			const auto& value = adPack[D].out_grad.val();  
			const size_t n_nonzero = adPack[D].out_grad.nnz();

			for (size_t k = 0; k < n_nonzero; ++k) {
				gradHere[index[k]] += value[k];
			}
#ifdef DEBUG
			Rcpp::Rcout << ",";
#endif			

		}
#ifdef DEBUG
		Rcpp::Rcout << "-\n";
#endif	
  #  pragma omp critical 
		{
			f += fOutHere;
			for(size_t Dpar=0;Dpar<n;Dpar++) {
				g[Dpar] += gradHere[Dpar];
			}
		}
	} // parallel
}




struct AD_Func_Opt {
	std::vector<GroupPack> &tape;
	Eigen::SparseMatrix<int> Htemplate;

    // ctor with no hessian pattern -> use empty storage
	AD_Func_Opt(std::vector<GroupPack> &tape_)
	: tape(tape_),
          Htemplate()     // bind ref to empty
          {}

    // ctor with hessian pattern
          AD_Func_Opt(std::vector<GroupPack> &tape_,
                std::vector<std::vector<int>> &hessianIJ) // i, p, x, lower triangle
          : tape(tape_), Htemplate()
          {
        // sanity check
          	if (hessianIJ.size() < 3) {
          		throw std::runtime_error("hessianIJ must have at least 3 components: row, colPtr, values");
          	}

        const int n    = static_cast<int>(hessianIJ[1].size()) - 1; // ncols (and nrows if square)
        const int nnz  = static_cast<int>(hessianIJ[0].size());

        // Map over the pattern data
        Eigen::Map<
        const Eigen::SparseMatrix<int, Eigen::ColMajor, int>
        > Hmap(
        	n, n, nnz,
            hessianIJ[1].data(),   // outerIndexPtr (p)
            hessianIJ[0].data(),   // innerIndexPtr (i)
            hessianIJ[2].data()    // valuePtr (x)
            );

        // copy pattern into H (so H owns its memory)
        Htemplate = Hmap;
    }

  template<class DerivedX>
    void get_hess(const Eigen::MatrixBase<DerivedX> &x,
    	Eigen::SparseMatrix<double> &H) {

    	using Index = Eigen::Index;

    	const Index n         = static_cast<Index>(x.size());
    	const Index NnonZero  = static_cast<Index>(Htemplate.nonZeros());
    	std::vector<double> outAll(NnonZero);

    	Eigen::Map<const Eigen::VectorXi> remapVec(
    		Htemplate.valuePtr(),
    		Htemplate.nonZeros()
    		);

    	CppAD::vector<double> xp(int(x.size()));
    	for (std::size_t i = 0; i < n; ++i) {
    		xp[i] = x[i];
    	}

 #pragma omp parallel 
    	{
    		CppAD::vector<double> w(1);
    		w[0] = 1.0;
    		const int Ngroups = tape.size();
    		std::vector<double> outHereAll(NnonZero, 0.0);

      #pragma omp for nowait
    		for (int g = 0; g < Ngroups; ++g) {
    			GroupPack &gp = tape[g];

    			gp.fun.sparse_hes(
    				xp,
                w,                         // weighting of outputs
                gp.out_hess,               // output container (sparse_rcv)
                gp.pattern_hess,           // sparsity pattern
                HESS_COLOR,                // typically "cppad.symmetric"
                gp.work_hess               // work space
                );
    			const std::vector<size_t>& matchHere = gp.match_hess;
    			const size_t Nhere = matchHere.size();
    			const CppAD::vector<double>& hessianOutHere = gp.out_hess.val();
    			for(size_t D=0;D < Nhere; D++) {
    				const size_t indexHere = matchHere[D];
    				outHereAll[indexHere] += hessianOutHere[D];
    			}
    		}

          #  pragma omp critical 
    		{	
    			for(size_t Dpar=0;Dpar<NnonZero;Dpar++) {
    				outAll[Dpar] += outHereAll[remapVec[Dpar]];
    			}
    		}
    } // parallel

    // 1. Make sure H has the same sparsity pattern as Htemplate
    if (H.rows()    != Htemplate.rows() ||
    	H.cols()    != Htemplate.cols() ||
    	H.nonZeros()!= Htemplate.nonZeros())
    {
        // Copy structure (and values) from template;
        // we'll overwrite the numeric values below.
    	H = Htemplate.cast<double>();
    }
    H.makeCompressed();

    Eigen::Map<Eigen::VectorXd>(H.valuePtr(), NnonZero) =
    Eigen::Map<const Eigen::VectorXd>(outAll.data(), NnonZero);
}

    // f and g together
    template <class DerivedX, class DerivedG>
void get_fdf(const Eigen::MatrixBase<DerivedX> &x,
	double &f,
	Eigen::MatrixBase<DerivedG> &g) {
	const std::size_t n = static_cast<std::size_t>(x.size());
	std::vector<double> gOut(n,0.0);
	CppAD::vector<double> xp(n);
	for (std::size_t i = 0; i < n; ++i) {
		xp[i] = x[i];
	}

	grad(xp, f, gOut, tape);
	for (std::size_t i = 0; i < n; ++i) {
		g[i] = gOut[i];
	}
}

    // f only
    template <class DerivedX>
void get_f(const Eigen::MatrixBase<DerivedX> &x,
	double &f) {

	const std::size_t n = static_cast<std::size_t>(x.size());
	CppAD::vector<double> xp(n);
	for (std::size_t i = 0; i < n; ++i) {
		xp[i] = x[i];
	}

	f = 0.0;
 #pragma omp parallel 
	{
		double fOutHere=0.0;

	      #pragma omp for nowait
		for(int D=0;D<tape.size();D++) {
			CppAD::vector<double> y = tape[D].fun.Forward(0, xp);
			fOutHere += y[0];
		}

  #  pragma omp critical 
		{
			f += fOutHere;
		}

	}
}

    // g only
    template <class DerivedX, class DerivedG>
void get_df(const Eigen::MatrixBase<DerivedX> &x,
	Eigen::MatrixBase<DerivedG> &g) {
	const std::size_t n = static_cast<std::size_t>(x.size());
	std::vector<double> gOut(n,0.0);
	CppAD::vector<double> xp(n);
	for (std::size_t i = 0; i < n; ++i) {
		xp[i] = x[i];
	}
	double f;
	grad(xp, f, gOut, tape);
	for (std::size_t i = 0; i < n; ++i) {
		g[i] = gOut[i];
	}
}
};
#endif