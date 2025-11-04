
#include"hpol.hpp"

const double sqrtDblEpsilon = DBL_EPSILON;//std::sqrt(DBL_EPSILON);


//#define DEBUG



// to do: put Q in adpack

double jointLogDens(
  std::vector<double> parameters, 
  Data& data, 
  Config& config
  ) {

  auto latent = unpack_params(parameters, data, config);

  std::vector<double> logLik(config.num_threads);
  std::vector<double> Qpart(config.num_threads);
  CppAD::vector<double> gammaScaled(data.Ngamma);

  const size_t NqDiag = data.Qdiag.size();

  omp_set_num_threads(config.num_threads);

      #pragma omp parallel
  { 

    const int tid=omp_get_thread_num();
    double logLikT = double(0.0), QpartT=double(0.0);

    #pragma omp for nowait
    for (size_t Dstrata = 0; Dstrata < data.Nstrata; Dstrata++) {
/*      loglikT += loglikOneStrata<double>(
        Dstrata,
        latent.gamma,
        latent,
        dat, 
        cfg
        )[0];*/

  auto etaHere = compute_eta_for_stratum<double, double>(
    Dstrata, data, latent.gamma, latent.beta);

  auto contrib = accumulate_contrib_for_stratum<double, double>(
    Dstrata, data, etaHere, latent, config
    );

    logLikT += contrib[0];

      }

  // Q diag.  
        #pragma omp for nowait
      for(size_t D=0;D<NqDiag;D++) {
        size_t mapHere = data.map[D];

        gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
        QpartT += latent.logTheta[mapHere] + gammaScaled[D]*gammaScaled[D]*(0.5*data.Qdiag[D]);
      }
      logLik[tid] = logLikT;      
      Qpart[tid] = QpartT;
  } // end parallel block


  double resultL=0.0, resultQ = 0.0;
  for(size_t D=0;D<config.num_threads;D++) {
//    Rcpp::Rcout << "com " << D << " " << logLik[D] << " " << Qpart[D] << "\n";
    resultL += logLik[D];
    resultQ += Qpart[D];
  }
  double result = -resultL + resultQ;

  // Q offdiag    
  for(size_t D = 0; D < data.Nq; D++) {
    result += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] * data.QsansDiag.x[D];
  }
  if (config.verbose ) Rcpp::Rcout << "L " << resultL << " Q " << resultQ << 
    " total " << result << "\n";

  return result;
}





Rcpp::NumericMatrix hessianQdense(
  const std::vector<double> parameters, 
  const Data& data,
  const Config& config) {


  const size_t Nparams = parameters.size();
  Rcpp::NumericMatrix hessianOut(Nparams, Nparams);
  auto fun = adFunQ(parameters, data, config);
  fun.Forward(0, parameters);

  const std::vector<double> w(1, 1.0);
  std::vector<double> u(Nparams, 0.0);

  for (int Dcol = 0; Dcol < Nparams; Dcol++) {
    std::fill(u.begin(), u.end(), 0.0);
    u[Dcol] = 1.0;
    fun.Forward(1, u);
    const auto ddw = fun.Reverse(2, w);
    for (int Drow = 0; Drow < Nparams; ++Drow) {
      hessianOut(Drow, Dcol) = ddw[2 * Drow + 1];
    }
  } // Dcol
  return(hessianOut);
}



std::vector<double> hessianQsparse(
      const std::vector<double> parameters, 
      const Data& data, 
      const Config& config
      ) {

      GroupPack Qpack = getAdFunQ(parameters, data, config);

      const size_t Nout = Qpack.outRowCol[0].size();
      std::vector<double> hessian(Nout);

      const std::vector<double> w(1, 1.0);

      Qpack.fun.SparseHessian(parameters, w, 
        Qpack.pattern, Qpack.outRowCol[0], 
        Qpack.outRowCol[1], hessian, Qpack.work);

      return(hessian);
}




Rcpp::S4 hessian(
  const std::vector<double>& parameters,
  std::vector<GroupPack>& adpack, 
  const Data& data,
  const Config& config
  ) {

  const Rcpp::List sparsity = config.group_sparsity;
  const size_t Nparams = parameters.size();
  const bool onlyRandom = Nparams == data.Ngamma;
  const size_t Ngroup = adpack.size();

  std::vector<double> qRes;
  std::vector<std::vector<double>> hessianOut(Ngroup);

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );


  #pragma omp parallel
  {

    const std::vector<double> w(1, 1.0);

  #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      const size_t Nhere = adpack[Dgroup].outRowCol[0].size();
      hessianOut[Dgroup] = std::vector<double>(Nhere);

//      adpack[Dgroup].fun.Forward(0, parameters);

      adpack[Dgroup].fun.SparseHessian(parameters, w, 
        adpack[Dgroup].pattern, adpack[Dgroup].outRowCol[0], 
        adpack[Dgroup].outRowCol[1], hessianOut[Dgroup], 
        adpack[Dgroup].work);
      } //Dgroup

  // add Q
  #pragma omp single 
    {
      // Q likelihood
      qRes = hessianQsparse(parameters, data, config);
    }

    } //parallel

// assemble
    Rcpp::S4 result=assembleHessian(hessianOut, qRes, sparsity, config, onlyRandom);


    return(result);
}


Rcpp::S4 hessian(
  const std::vector<double>& parameters,
  const Data& data,
  const Config& config
  ) {

  const Rcpp::List strata = config.groups;
  const Rcpp::List sparsity = config.group_sparsity;

  std::vector<GroupPack> adpack = getAdFun(
    parameters, data, config);

  Rcpp::S4 result = hessian(
    parameters, adpack, data, config);
  return(result);
}


  Rcpp::LogicalMatrix hessianDenseLogical(
    const std::vector<double> parameters,
    const Data& data,
    const Config& config
    ) {

  const Rcpp::List strata = config.groups;
  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

  const size_t Nparams = parameters.size();
  const size_t Ngroup = strataP.size()-1;

  Rcpp::LogicalMatrix result(Nparams, Nparams);


  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {

    const std::vector<double> w(1, 1.0);
    std::vector<double> u(Nparams, 0.0);

    std::vector<std::vector<bool>> hessianOutHere(Nparams);
    for(size_t D=0;D<Nparams;++D){
      hessianOutHere[D] = std::vector<bool>(Nparams, false);
    }

  #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      auto fun=adFunGroup(
        parameters, data, config, strataI, 
        strataP[Dgroup], strataP[Dgroup+1]);
      fun.Forward(0, parameters);

      for (int Dcol = 0; Dcol < Nparams; Dcol++) {
        std::vector<bool>& hessianOutHereCol = hessianOutHere[Dcol];
        std::fill(u.begin(), u.end(), 0.0);
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        auto outHere = fun.Reverse(2, w);
        for (int Drow = 0; Drow < Nparams; ++Drow) {
          hessianOutHereCol[Drow] = hessianOutHereCol[Drow] || ( std::fabs(outHere[2 * Drow + 1]) > sqrtDblEpsilon) ;
        }
      } // Dcol
  } // Dgoup

 #pragma omp single 
  {
      // Q likelihood
    const bool useQ = data.Qdiag.size();
    if(useQ) {
      auto fun = adFunQ(parameters, data, config);
      fun.Forward(0, parameters);
      const std::vector<double> w(1, 1.0);
      std::vector<double> u(Nparams, 0.0);

      for (int Dcol = 0; Dcol < Nparams; Dcol++) {
        std::fill(u.begin(), u.end(), 0.0);
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        const auto ddw = fun.Reverse(2, w);
        std::vector<bool>& hessianOutHereCol = hessianOutHere[Dcol];

        for (int Drow = 0; Drow < Nparams; ++Drow) {
          hessianOutHereCol[Drow] = hessianOutHereCol[Drow] || ( fabs(ddw[2 * Drow + 1]) > sqrtDblEpsilon);
        }
      } // Dcol
    } // useQ
  } // single

#pragma omp critical
{
for(size_t Dcol=0;Dcol<Nparams;Dcol++) {
  std::vector<bool>& hessianOutHereCol = hessianOutHere[Dcol];
  for(size_t Drow=0;Drow<Nparams;Drow++) {
    result(Drow, Dcol) = result(Drow, Dcol) || hessianOutHereCol[Drow];
  }
}  
} // critical
} // parallel

return(result);
}

  Rcpp::NumericMatrix hessianDense(
    const std::vector<double>& parameters,
    const Data& data,
    const Config& config
    ) {

  const Rcpp::List strata=config.groups;
  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

  const size_t Nparams = parameters.size();
  const size_t Ngroup = strataP.size()-1;
  const bool useQ = data.Qdiag.size()>0;


  Rcpp::NumericMatrix result(Nparams, Nparams);


  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
    CppAD::thread_alloc::hold_memory(true);

    const std::vector<double> w(1, 1.0);
    std::vector<double> u(Nparams, 0.0);

    std::vector<std::vector<double>> hessianOutHere(Nparams);
    for(size_t D=0;D<Nparams;++D){
      hessianOutHere[D] = std::vector<double>(Nparams, 0.0);
    }

  #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      auto fun=adFunGroup(
        parameters, data, config, strataI, 
        strataP[Dgroup], strataP[Dgroup+1]);
      fun.Forward(0, parameters);

      for (size_t Dcol = 0; Dcol < Nparams; Dcol++) {
        std::vector<double>& hessianOutHereCol = hessianOutHere[Dcol];
        std::fill(u.begin(), u.end(), 0.0);
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        auto outHere = fun.Reverse(2, w);
        for (int Drow = 0; Drow < Nparams; ++Drow) {
          hessianOutHereCol[Drow] += outHere[2 * Drow + 1];
        }
      } // Dcol
  } // Dgoup

 #pragma omp single 
  {
      // Q likelihood
    if(useQ) {
      auto fun = adFunQ(parameters, data, config);
      fun.Forward(0, parameters);
      const std::vector<double> w(1, 1.0);
      std::vector<double> u(Nparams, 0.0);

      for (size_t Dcol = 0; Dcol < Nparams; Dcol++) {
        std::fill(u.begin(), u.end(), 0.0);
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        const auto ddw = fun.Reverse(2, w);

        std::vector<double>& hessianOutHereCol = hessianOutHere[Dcol];

        for (size_t Drow = 0; Drow < Nparams; ++Drow) {
          hessianOutHereCol[Drow] += ddw[2 * Drow + 1];
        }
      } // Dcol
    } // useQ
  } // single

#pragma omp critical
{
for(size_t Dcol=0;Dcol<Nparams;Dcol++) {
  std::vector<double>& hessianOutHereCol = hessianOutHere[Dcol];
  for(size_t Drow=0;Drow<Nparams;Drow++) {
    result(Drow, Dcol) +=  hessianOutHereCol[Drow];
  }
}  
} // critical
} // parallel

return(result);
}




std::vector<double> grad(
  const std::vector<double> parameters,
  std::vector<GroupPack>& adpack, 
  const Data& data,
  const Config& config
  ) {

  const size_t Nparams = parameters.size();
  const size_t Ngroup = adpack.size();
  const bool useQ = data.Qdiag.size()>0;

  std::vector<double> gradOut(Nparams, 0);

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
    CppAD::thread_alloc::hold_memory(true); 
    std::vector<double> gradHere(Nparams, 0);
    std::vector<double> w(1, 1.0);

      #pragma omp for nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      adpack[Dgroup].fun.Forward(0, parameters);
      auto gradThisGroup = adpack[Dgroup].fun.Reverse(1, w);

      for(size_t D=0;D<Nparams;D++) {
        gradHere[D] += gradThisGroup[D];
      }

    } // group

// get rid?
#pragma omp barrier

#pragma omp single 
    {
      // Q likelihood
      if(useQ) {
        auto fun = adFunQ(parameters, data, config);
        fun.Forward(0, parameters);
        auto gradQ = fun.Reverse(1, w);
        for(size_t D=0;D<Nparams;D++) {
          gradHere[D]+= gradQ[D];
        }
      } // useQ
    } // single Q

#pragma omp critical
{
for(size_t D=0;D<Nparams;D++) {
  gradOut[D] += gradHere[D];
}  

} // critical
  }// parallel


  return(gradOut);

}

std::vector<double> gradNeedAd(
  const std::vector<double>& parameters,
  const Data& data,
  const Config& config
  ) {

  std::vector<GroupPack> adpack = getAdFun(
    parameters, data, config);

  std::vector<double>  result = grad(
    parameters, adpack, data, config);
  return(result);
}




/* R exported stuff */
//' @export
// [[Rcpp::export]]
SEXP getAdFun(
  Rcpp::NumericVector x, 
  const Rcpp::List data, 
  const Rcpp::List config
  ) {

  Data   dataC(data);
  Config configC(config);
   std::vector<double> parametersC = Rcpp::as<std::vector<double>>(x);


  std::vector<GroupPack> adpack = getAdFun(
    parametersC, dataC, configC);

  auto* ptr = new std::vector<GroupPack>(std::move(adpack));

  Rcpp::XPtr<std::vector<GroupPack>> xp(ptr, /*deleteOnFinalizer=*/true);

  // (Optional) tag a class so you can validate on the R side
  xp.attr("class") = "adpack_ptr";
  return xp;
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
 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);


  double result= jointLogDens(
    parametersC,
    dat, cfg);

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

  std::vector<double> parametersC(parameters.begin(), parameters.end());
  const Data   dataC(data);
  const Config configC(config);

  AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
  std::vector<GroupPack>* packsPtr = ad.ptr;

  std::vector<double> result = grad(parametersC, *packsPtr, dataC, configC);
  ad.cleanup();

  return Rcpp::NumericVector(result.begin(), result.end());
}

//' @export
// [[Rcpp::export]]
Rcpp::S4 hessian(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config,
  SEXP adFun = R_NilValue)
{

 const Data dataC(data);
 const Config configC(config);

 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);

AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
std::vector<GroupPack>* packsPtr = ad.ptr;


 Rcpp::S4 result = hessian(
  parametersC, *packsPtr, dataC, configC
  );
  ad.cleanup();

 return(result);
}


//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hessianQdense(
      const Rcpp::NumericVector parameters, 
      const Rcpp::List& data,
      const Rcpp::List &config) {

      const Data dataC(data);
      const Config configC(config);

      std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);
      Rcpp::NumericMatrix result = hessianQdense(parametersC, dataC, configC);
      return(result);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hessianDense(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config
  ) {

 const Data dataC(data);
 const Config configC(config);

 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);

Rcpp::NumericMatrix result = hessianDense(
  parametersC, dataC, configC
  );
 return(result);
}

//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix hessianDenseLogical(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config
  ) {

 const Data dataC(data);
 const Config configC(config);

 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);

Rcpp::LogicalMatrix result = hessianDenseLogical(
  parametersC, dataC, configC
  );
 return(result);
}



//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix gradLogical(
 const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config
  ) {

    const Data dataC(data);
    const Config configC(config);

    std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);
    const size_t Nparams = parametersC.size();

Rcpp::IntegerVector idx(Nparams);
std::iota(idx.begin(), idx.end(), 0);

    Rcpp::LogicalMatrix result(Nparams, dataC.Nstrata);
         std::vector<double> w(1, 1.0);

  omp_set_num_threads(configC.num_threads);
  CppAD::thread_alloc::parallel_setup(
    configC.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
      CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
      CppAD::AD<double> loglik = 0;
      for (size_t D = 0; D < Nparams; D++) {
        ad_params[D] = parametersC[D];  
      }

      #pragma omp for nowait
    for(size_t Dstrata = 0;Dstrata< dataC.Nstrata;Dstrata++) {


      CppAD::Independent(ad_params);

      auto latent=unpack_params(ad_params, dataC, configC);


        auto etaHere = compute_eta_for_stratum(
          Dstrata, dataC, latent.gamma, latent.beta);

        auto out = accumulate_contrib_for_stratum(
          Dstrata, dataC, etaHere, latent, configC);

      CppAD::ADFun<double> fun(ad_params, out);

      fun.Forward(0, parametersC);
      auto gradHere = fun.Reverse(1, w);
      for(size_t D=0;D<Nparams;D++) {
        result(D, Dstrata) = (fabs(gradHere[D])>sqrtDblEpsilon);
      }  
    }
  } // parallel
    return(result);
}
