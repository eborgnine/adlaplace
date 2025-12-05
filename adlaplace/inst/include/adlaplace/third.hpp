#ifndef THIRD_HPP
#define THIRD_HPP

inline CppAD::vector<double> traceHinvT(
    const CppAD::vector<double>&  parameters,
    const CSCMatrix& LinvPt,
    const CSCMatrix& LinvPtColumns,
    const Config& config,
    std::vector<GroupPack>& fun
    ) {

    const size_t Nparams = parameters.size();
    const size_t NcolsL = LinvPt.ncol;
    const size_t Nbeta = config.beta.size();
    const size_t Ngroup = fun.size();
    CppAD::vector<double> result(Nparams, 0.0);


  #pragma omp parallel
    {
        CppAD::vector<double> w{0.0, 0.0, 1.0};  
        CppAD::vector<double> x_val = parameters; 
        CppAD::vector<double> directionZeros(Nparams, 0.0);
        CppAD::vector<double> direction(Nparams, 0.0);
        CppAD::vector<double> outHere(Nparams, 0.0);

    # pragma omp for nowait
        for(int Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
            const int colEnd = LinvPtColumns.p[Dgroup+1];

            if(config.verbose) {
                Rcpp::Rcout << " " << Dgroup << " " << colEnd  << " \n";
            }        

            for( int Dp = LinvPtColumns.p[Dgroup]; Dp < colEnd; Dp++) {
              const int Dcol = LinvPtColumns.i[Dp];
              std::fill(direction.begin(), direction.end(), 0.0);
              const int endHere = LinvPt.p[Dcol+1];
              for(int D = LinvPt.p[Dcol];D<endHere;D++) {
                const int Ihere = LinvPt.i[D];
                direction[Ihere + Nbeta] = LinvPt.x[D];            
            }
            fun[Dgroup].fun.Forward(0, x_val);
            fun[Dgroup].fun.Forward(1, direction);
            fun[Dgroup].fun.Forward(2, directionZeros);
            auto dw = fun[Dgroup].fun.Reverse(3, w);
            for(size_t Dk=0;Dk<Nparams;Dk++) {
                outHere[Dk] += dw[3*Dk];
            } 
        } // Dgamma
    }   // Dgroup

        #pragma omp critical(accum) 
    {
        for (int Dk = 0; Dk < Nparams; ++Dk) {
            result[Dk] += outHere[Dk];
        }
        if(config.verbose) {
            Rcpp::Rcout << " " << omp_get_thread_num()  << " \n";
        }        

    } // critical
    }// parallel
    return(result);
}

#endif

