
#include"hpol.hpp"

//#define DEBUG


/*#if defined(CPPAD_HAS_COLPACK) && CPPAD_HAS_COLPACK
static const char* JAC_COLOR = "colpack";
static const char* HESS_COLOR = "colpack.symmetric";
#else*/
  static const char* JAC_COLOR = "cppad";  // fallback if ColPack not available
  static const char* HESS_COLOR = "cppad.symmetric";
//#endif



double jointLogDens( //is this slower?  woudl be convenient if not
  const CppAD::vector<double> parameters,
  std::vector<GroupPack>& adpack
  ) {
  const size_t Ngroup = adpack.size();
  double result=0.0;

#ifdef DEBUG
  Rcpp::Rcout << "here ";
#endif
  #pragma omp parallel
  {

    CppAD::vector<double> x = parameters;
    double resultHere = 0.0;

    #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
#ifdef DEBUG
  Rcpp::Rcout << "G " << Dgroup << " ";
#endif

      auto resultHere1 = adpack[Dgroup].fun.Forward(0, x);
      resultHere += resultHere1[0];
    }
      #  pragma omp critical 
    {
      result += resultHere;
    }

  }
#ifdef DEBUG
  Rcpp::Rcout << "here2 \n";
#endif

  return(result);
}

void grad(
  const CppAD::vector<double> parameters,
  std::vector<GroupPack>& adpack, 
  std::vector<double>& result
  ) {

  const size_t Nparams = parameters.size();
  const size_t Ngroup = adpack.size();

  #pragma omp parallel
  {

    CppAD::vector<double> x = parameters;
    std::vector<double> gradHere(Nparams, 0);
    CppAD::vector<double> w(1);  
    w[0] = 1.0;

      #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
      auto& gp = adpack[Dgroup];

//      gp.fun.Forward(0,x);
      gp.fun.sparse_jac_rev(
        x,
        gp.out_grad,         // CppAD::sparse_rcv<...>, values overwritten
        gp.pattern_grad,
        JAC_COLOR,      
        gp.work_grad);

    const auto& index = gp.out_grad.col();   // indices of nonzero gradient components
    const auto& value = gp.out_grad.val();   // corresponding values
    const size_t n_nonzero = gp.out_grad.nnz();

    for (size_t k = 0; k < n_nonzero; ++k) {
      gradHere[index[k]] += value[k];
    }
    } // group

  #  pragma omp critical 
  // pragma parallel reduction(+: out[:n])
    {
      for (int d = 0; d < Nparams; ++d) {
        result[d] += gradHere[d];
      }
    }

  }// parallel

}


void hessian(
  const CppAD::vector<double>& parameters,
  std::vector<GroupPack>& adpack, 
  std::vector<double>& result
  ) {

  const size_t Ngroup = adpack.size(); // includes Q
  const size_t NnonZero = result.size();

  #pragma omp parallel
  {

    CppAD::vector<double> w(1);
    w[0] = 1.0;
    std::vector<double> outHereAll(NnonZero, 0.0);

  #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      adpack[Dgroup].fun.sparse_hes(parameters, w, 
        adpack[Dgroup].out_hess, adpack[Dgroup].pattern_hess, 
        HESS_COLOR, adpack[Dgroup].work_hess);

      const std::vector<size_t>& matchHere = adpack[Dgroup].outRowCol[2];
      const size_t Nhere = matchHere.size();
      const CppAD::vector<double>& hessianOutHere = adpack[Dgroup].out_hess.val();

      for(size_t D=0;D < Nhere; D++) {
        const size_t indexHere = matchHere[D];
        outHereAll[indexHere] += hessianOutHere[D];
      }
    } //Dgroup

#pragma omp critical
    {
      for(size_t D=0;D<NnonZero;D++) {
        result[D] +=  outHereAll[D];
      }  
    } // critical

  } //parallel
}

// converts result into a gCMatrix
Rcpp::RObject hessian(
  const CppAD::vector<double>& parameters,
  std::vector<GroupPack>& adpack, 
  const Rcpp::List sparsity
  ) {

  const Rcpp::IntegerVector pOut = sparsity["p"];
  const Rcpp::IntegerVector iOut = sparsity["i"];
  const Rcpp::IntegerVector jOut = sparsity["j"];
  const int Nnonzero = iOut.size();
  const int N = adpack[0].out_hess.nr();

  std::vector<double> hessianOut(Nnonzero);

  hessian(parameters, adpack, hessianOut);
  auto result = make_convert_gCmatrix(hessianOut, iOut, jOut, N);

  return(result);
}


// work out random or full sizes in adpack, pass approprite sparsity
Rcpp::RObject hessian(
  const CppAD::vector<double>& parameters,
  std::vector<GroupPack>& adpack, 
  const Config& config
  ) {

  const int N = adpack[0].out_hess.nr();

  const Rcpp::List secondAll = config.sparsity["second"];
  const Rcpp::List secondRandom = secondAll["random"];
  const Rcpp::List secondFull = secondAll["full"];
  const Rcpp::IntegerVector pRandom = secondRandom["p"];
  const Rcpp::IntegerVector pFull= secondFull["p"];
  const int Nrandom = pRandom.size()-1, Nfull = pFull.size()-1;
  const bool onlyRandom = Nrandom == N;

  const Rcpp::List sparsityOut = onlyRandom?secondRandom:secondFull;

  auto result = hessian(parameters, adpack, sparsityOut);

  return(result);
}



//' @export
// [[Rcpp::export]]
double jointLogDens(
  Rcpp::NumericVector parameters, 
  const Rcpp::List data, 
  const Rcpp::List config,
  SEXP adFun = R_NilValue)
{

  Data   dat(data);
  Config cfg(config);

  const size_t N = parameters.size();
  CppAD::vector<double> parametersC(N);
  for(size_t D=0; D<N;D++) {
    parametersC[D] = parameters[D];
  }

  omp_set_num_threads(cfg.num_threads);
  CppAD::thread_alloc::parallel_setup(
    cfg.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  AdpackHandle ad = getAdpackFromR(adFun, parametersC, dat, cfg);
  std::vector<GroupPack>* packsPtr = ad.ptr;

  double result= jointLogDens(parametersC, *packsPtr);

  return(result);

}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grad(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config,
  SEXP adFun = R_NilValue)
{

  const size_t N = parameters.size();
  CppAD::vector<double> parametersC(N);
  for(size_t D=0; D<N;D++) {
    parametersC[D] = parameters[D];
  }
  const Data   dataC(data);
  const Config configC(config);

  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
  std::vector<GroupPack>* packsPtr = ad.ptr;
  std::vector<double> result(N);

  grad(parametersC, *packsPtr, result);

  ad.cleanup();

  return Rcpp::NumericVector(result.begin(), result.end());
}



//' @export
// [[Rcpp::export]]
Rcpp::RObject hessian(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config,
  SEXP adFun = R_NilValue)
{

 const Data dataC(data);
 const Config configC(config);

 const size_t N = parameters.size();
 CppAD::vector<double> parametersC(N);
 for(size_t D=0; D<N;D++) {
  parametersC[D] = parameters[D];
}

  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
std::vector<GroupPack>* packsPtr = ad.ptr;

Rcpp::RObject result = hessian(
  parametersC, *packsPtr, configC
  );

ad.cleanup();

return(result);
}


