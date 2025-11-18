#include"hpol.hpp"


// Simple CSC container suitable for OpenMP use
struct CSCMatrix {
    int nrow;
    int ncol;
    std::vector<int>    p;  // column pointers (size ncol + 1)
    std::vector<int>    i;  // row indices (size nnz)
    std::vector<double> x;  // values      (size nnz)
};

// Construct CSCMatrix from an R "dgCMatrix"
inline CSCMatrix makeCSC(const Rcpp::S4& M)
{
    CSCMatrix out;

    Rcpp::IntegerVector Dim = M.slot("Dim");
    Rcpp::IntegerVector pR  = M.slot("p");
    Rcpp::IntegerVector iR  = M.slot("i");
    Rcpp::NumericVector xR  = M.slot("x");

    out.nrow = Dim[0];
    out.ncol = Dim[1];

    const int nnz = xR.size();

    out.p.assign(pR.begin(), pR.end());
    out.i.assign(iR.begin(), iR.end());
    out.x.assign(xR.begin(), xR.end());

    return out;
}


CppAD::vector<double> traceHinvT(
    const CppAD::vector<double>&  parameters,
    const CSCMatrix& LinvPt,
    const Config& config,
    std::vector<GroupPack>& fun,
    GroupPack& Qfun
    ) {

    const size_t Nparams = parameters.size();
    const size_t NcolsL = LinvPt.ncol;
    const size_t Nbeta = config.beta.size();
    const size_t Ngroup = fun.size();
    CppAD::vector<double> result(Nparams, 0.0);

    omp_set_num_threads(config.num_threads);
    CppAD::thread_alloc::parallel_setup(
        config.num_threads,
        [](){ return in_parallel_wrapper(); },
        [](){ return static_cast<size_t>(thread_num_wrapper()); }
        );

  #pragma omp parallel
    {
        CppAD::vector<double> w{0.0, 0.0, 1.0};  
        CppAD::vector<double> x_val = parameters; 
        CppAD::vector<double> directionZeros(Nparams, 0.0);
        CppAD::vector<double> direction(Nparams, 0.0);
        CppAD::vector<double> outHere(Nparams, 0.0);

    # pragma omp for nowait
        for(int Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
            for( int Dcol = 0; Dcol < NcolsL; Dcol++) {

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
    # pragma omp single
    {
        for( int Dcol = 0; Dcol < NcolsL; Dcol++) {

          std::fill(direction.begin(), direction.end(), 0.0);
          const int endHere = LinvPt.p[Dcol+1];
          for(int D = LinvPt.p[Dcol];D<endHere;D++) {
            const int Ihere = LinvPt.i[D];
            direction[Ihere + Nbeta] = LinvPt.x[D];            
          }
        Qfun.fun.Forward(0, x_val);
        Qfun.fun.Forward(1, direction);
        Qfun.fun.Forward(2, directionZeros);
        auto dw = Qfun.fun.Reverse(3, w);
        for(size_t Dk=0;Dk<Nparams;Dk++) {
            outHere[Dk] += dw[3*Dk];
        } 
        } // Dgamma

    } // single

        #pragma omp critical(accum)
    for (int Dk = 0; Dk < Nparams; ++Dk) {
        result[Dk] += outHere[Dk];
    }

    }// parallel
    return(result);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector traceHinvT(
  const Rcpp::NumericVector parameters,
  const Rcpp::S4& LinvPt,
  const Rcpp::List data, 
  const Rcpp::List config,
  SEXP adFun = R_NilValue
  ) {

 const Data   dataC(data);
 const Config configC(config);
 const size_t Nparams = parameters.size();
 CppAD::vector<double> parametersC(Nparams);

 for(size_t D=0; D<Nparams; D++) {
    parametersC[D] = parameters[D];
}

CSCMatrix LinvPtC = makeCSC(LinvPt);

AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
std::vector<GroupPack>* fun = ad.ptr;
auto Qfun = getAdFunQ(parametersC, dataC, configC);

auto resultC = traceHinvT(parametersC, LinvPtC, configC, *fun, Qfun);

Rcpp::NumericVector result(Nparams);
for(size_t D=0; D<Nparams; D++) {
    result[D] = resultC[D];
}
return(result);
}

CppAD::vector<double> third_quadratic(
    const CppAD::vector<double>&  parameters, // beta, gamma, theta
    const CppAD::vector<double>& direction, // x
    const Config& config,
    std::vector<GroupPack>& fun,
    GroupPack& Qfun
    ) {
    const size_t Nparams = parameters.size();
    const size_t Ngroup = fun.size();
    CppAD::vector<double> result(Nparams, 0.0);

    omp_set_num_threads(config.num_threads);
    CppAD::thread_alloc::parallel_setup(
        config.num_threads,
        [](){ return in_parallel_wrapper(); },
        [](){ return static_cast<size_t>(thread_num_wrapper()); }
        );

  #pragma omp parallel
    {
        const CppAD::vector<double> w{0.0, 0.0, 1.0};  
        CppAD::vector<double> x_val = parameters; 
        const CppAD::vector<double> directionZeros(Nparams, 0.0);
        CppAD::vector<double> outHere(Nparams, 0.0);

    # pragma omp for nowait
        for(int Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
          fun[Dgroup].fun.Forward(0, x_val);
          fun[Dgroup].fun.Forward(1, direction);
          fun[Dgroup].fun.Forward(2, directionZeros);
          auto dw = fun[Dgroup].fun.Reverse(3, w);
          for(size_t Dk=0;Dk<Nparams;Dk++) {
            outHere[Dk] += dw[3*Dk];
        } 
    }
    # pragma omp single
    {
        Qfun.fun.Forward(0, x_val);
        Qfun.fun.Forward(1, direction);
        Qfun.fun.Forward(2, directionZeros);
        auto dw = Qfun.fun.Reverse(3, w);
        for(size_t Dk=0;Dk<Nparams;Dk++) {
            outHere[Dk] += dw[3*Dk];
        } 
    } // single

    #pragma omp critical(accum)
    for (int Dk = 0; Dk < Nparams; ++Dk) {
        result[Dk] += outHere[Dk];
    }

  } // parallel

  return(result);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector thirdPolarization(
  const Rcpp::NumericVector parameters,
  const Rcpp::NumericVector A, 
  const Rcpp::NumericVector B, 
  const Rcpp::List data, 
  const Rcpp::List config,
  SEXP adFun = R_NilValue
  ) {


 const Data   dataC(data);
 const Config configC(config);
 const size_t Nparams = parameters.size();
 CppAD::vector<double> parametersC(Nparams), ApB(Nparams), AmB(Nparams);


 for(size_t D=0; D<Nparams; D++) {
    parametersC[D] = parameters[D];
    ApB[D] = A[D] + B[D];
    AmB[D] = A[D] - B[D];
}

AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
std::vector<GroupPack>* fun = ad.ptr;
auto Qfun = getAdFunQ(parametersC, dataC, configC);


// Q_plus[k] = sum_{ij} T_{ijk} (A_i + B_i)(A_j + B_j)
auto Qplus = third_quadratic(parametersC, ApB, configC, *fun, Qfun);

// Q_minus[k] = sum_{ij} T_{ijk} (A_i - B_i)(A_j - B_j)
auto Qminus = third_quadratic(parametersC, AmB, configC, *fun, Qfun);

// mixed[k] = sum_{ij} T_{ijk} A_i B_j
Rcpp::NumericVector result(Nparams);
for(size_t D=0; D<Nparams; D++) {
    result[D] = 0.5 * (Qplus[D] - Qminus[D]);
}

return(result);
}