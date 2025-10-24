
#include"hpol.hpp"
#include"loglikHelpers.hpp"
#include"matrixUtils.hpp"




//#define DEBUG

#include<omp.h>



//' @export
// [[Rcpp::export]]
double jointLogDens(
  std::vector<double> x, 
  Data& data, 
  Config& config
  ) {

  auto latent = unpack_params(x, data, config);

  std::vector<double> logLik(config.num_threads);
  std::vector<double> Qpart(config.num_threads);
  CppAD::vector<double> gammaScaled(data.Ngamma);


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
      for(size_t D=0;D<data.Ngamma;D++) {
        size_t mapHere = data.map[D];

        gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
        QpartT += latent.logTheta[mapHere] + gammaScaled[D]*gammaScaled[D]*(0.5*data.Qdiag[D]);
      }
      logLik[tid] = logLikT;      
      Qpart[tid] = QpartT;
  } // end parallel block


  double resultL=0.0, resultQ = 0.0;
  for(size_t D=0;D<config.num_threads;D++) {
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




/* likelihood part */
CppAD::ADFun<double> adFunGroup(
      const std::vector<double> & parameters,  
      const Data& data, 
      const Config& config,
      const Rcpp::IntegerVector& strataI,
      const size_t start,
      const size_t end
      ) {

      const size_t Nparams = parameters.size();
      CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
      CppAD::vector<CppAD::AD<double>> minusLogDens(1);
      CppAD::AD<double> loglik = 0;

      for (size_t D = 0; D < Nparams; D++) {
        ad_params[D] = parameters[D];  
      }
      CppAD::Independent(ad_params);

      auto latent=unpack_params(parameters, data, config);

      for (size_t Dindex = start; Dindex < end;  Dindex++) {

        const size_t Dstrata = strataI[Dindex];

        auto etaHere = compute_eta_for_stratum(
          Dstrata, data, latent.gamma, latent.beta);

        auto contrib = accumulate_contrib_for_stratum(
          Dstrata, data, etaHere, latent, config);

        loglik += contrib[0];
  } // Dstrata

  minusLogDens[0] = -loglik;

  CppAD::ADFun<double> fun(ad_params, minusLogDens);

  return(fun);
}




CppAD::ADFun<double> adFunQ(
  const std::vector<double> & parameters,  
  const Data& data,
  const Config& config) {

  const size_t Nparams = parameters.size();
  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
  for(size_t D=0;D<Nparams;D++) {
    ad_params[D] = parameters[D];
  }
  CppAD::Independent(ad_params);   


  auto latent=unpack_params(ad_params, data, config);

  CppAD::vector<CppAD::AD<double>> gammaScaled(data.Ngamma);
  CppAD::vector<CppAD::AD<double>> result(1);

  CppAD::AD<double> result0=0;
  const Rcpp::IntegerVector map = data.map;
  for (size_t D = 0; D < data.Ngamma; ++D) {
    const size_t mapHere = map[D];

    auto thetaHere = latent.theta[mapHere];
    auto logThetaHere = latent.logTheta[mapHere];

    gammaScaled[D] = latent.gamma[D] / thetaHere;
    result0 += logThetaHere + (0.5 * data.Qdiag[D]) * gammaScaled[D] * gammaScaled[D];
  }

      // Q offdiag    
  for(size_t D = 0; D < data.Nq; D++) {
    result0 += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] 
    * data.QsansDiag.x[D];
  }

  result[0] = result0;

  CppAD::ADFun<double> fun(ad_params, result);

  return(fun);
}




Rcpp::NumericMatrix hessianQdense(
  const std::vector<double> parameters, 
  const Data& data,
  const Config& config) {

  const size_t Nparams = parameters.size();
  Rcpp::NumericMatrix hessianOut(Nparams, Nparams);
  auto fun = adFunQ(parameters, data, config);
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
      const Config& config,
      const bool onlyRandom
      ) {

      const size_t Nparams = parameters.size();
      const Rcpp::List sparsity=config.sparsity;


      const Rcpp::List sparsityQ = sparsity["Q"]; 
      const Rcpp::List nsQ =  sparsityQ["nonSymmetric"]; 
      const Rcpp::List outList = onlyRandom?sparsityQ["random"]:sparsityQ["full"];
      

      const Rcpp::IntegerVector rowQns = nsQ["i"]; 
      const Rcpp::IntegerVector colQns = nsQ["j"]; 
      const auto pattern = build_pattern_from_R(rowQns, colQns, Nparams);

      const Rcpp::IntegerVector rowOut = outList["i"]; 
      const Rcpp::IntegerVector colOut = outList["j"]; 

      const std::vector<size_t> SrowOut = Rcpp::as<std::vector<size_t>>(rowOut);
      const std::vector<size_t> ScolOut = Rcpp::as<std::vector<size_t>>(colOut);
      const size_t Nout =  SrowOut.size();


      std::vector<double> hessian(Nout);

      CppAD::sparse_hessian_work work;
      const std::vector<double> w(1, 1.0);

      auto fun = adFunQ(parameters, data, config);
      fun.SparseHessian(parameters, w, pattern, SrowOut, ScolOut, hessian, work);

      return(hessian);
}






Rcpp::S4 hessian(
  const std::vector<double> parameters,
  const Data& data,
  const Config& config,
  const Rcpp::List& strata,
  const Rcpp::List& sparsity,
  const bool onlyRandom
  ) {

  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

  const size_t Nparams = parameters.size();
  const size_t Ngroup = strataP.size()-1;

  std::vector<std::vector<double>> hessianOut(Ngroup);
  std::vector<double> qRes;

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {

    const std::vector<double> w(1, 1.0);

  #pragma omp for
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      const Rcpp::List sparsityHere=sparsity[Dgroup];
      const Rcpp::List second = sparsityHere["second"];

        // pattern
      const Rcpp::List nonSymmetric = second["nonSymmetric"];
      const Rcpp::IntegerVector Srow = nonSymmetric["i"]; 
      const Rcpp::IntegerVector Scol = nonSymmetric["j"]; 
      auto pattern = build_pattern_from_R(Srow, Scol, Nparams);
      CppAD::sparse_hessian_work work;

        // out
      const Rcpp::List outList = onlyRandom?second["random"]:second["full"];
      const Rcpp::IntegerVector SrowOutR = outList["i"]; 
      const Rcpp::IntegerVector ScolOutR = outList["j"]; 
      const std::vector<size_t> SrowOut = Rcpp::as<std::vector<size_t>>(SrowOutR);
      const std::vector<size_t> ScolOut = Rcpp::as<std::vector<size_t>>(ScolOutR);

      auto fun=adFunGroup(
        parameters, data, config, strataI, 
        strataP[Dgroup], strataP[Dgroup+1]);

      std::vector<double> out(SrowOut.size());
      fun.SparseHessian(parameters, w, pattern, SrowOut, ScolOut, out, work);

      hessianOut[Dgroup] = out;
      } //Dgroup


  // add Q
  #pragma omp single 
    {
      // Q likelihood
      qRes = hessianQsparse(parameters, data, config, onlyRandom);
    }

    } //parallel

// assemble

    Rcpp::S4 result = assembleHessian(hessianOut, qRes, sparsity, config, onlyRandom);

    return(result);
}

  Rcpp::NumericMatrix hessianDense(
    const std::vector<double> parameters,
    const Data& data,
    const Config& config,
    const Rcpp::List& strata
    ) {


  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

  const size_t Nparams = parameters.size();
  const size_t Ngroup = strataP.size()-1;

  std::vector<std::vector<std::vector<double>>> hessianOut(config.num_threads);


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

    std::vector<std::vector<double>> hessianOutHere(Nparams);
    for(size_t D=0;D<Nparams;++D){
      hessianOutHere[D] = std::vector<double>(Nparams, 0);
    }

  #pragma omp for
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      auto fun=adFunGroup(
        parameters, data, config, strataI, 
        strataP[Dgroup], strataP[Dgroup+1]);

      fun.Forward(0, parameters);
      for (int Dcol = 0; Dcol < Nparams; Dcol++) {
        std::vector<double> hessianOutHereCol = hessianOutHere[Dcol];
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
    const bool useQ = data.Qdiag.size();
    if(useQ) {
      auto fun = adFunQ(parameters, data, config);
      const std::vector<double> w(1, 1.0);
      std::vector<double> u(Nparams, 0.0);

      for (int Dcol = 0; Dcol < Nparams; Dcol++) {
        std::fill(u.begin(), u.end(), 0.0);
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        const auto ddw = fun.Reverse(2, w);
        std::vector<double> hessianOutHereCol = hessianOutHere[Dcol];

        for (int Drow = 0; Drow < Nparams; ++Drow) {
          hessianOutHereCol[Drow] = ddw[2 * Drow + 1];
        }
      } // Dcol
    } // useQ
  } // single
  const int tid=omp_get_thread_num();
  hessianOut[tid] = hessianOutHere;
} // parallel

Rcpp::NumericMatrix result(Nparams, Nparams);
for(size_t Dthread=0;Dthread<config.num_threads;++Dthread) {
  std::vector<std::vector<double>> hessianOutHere = hessianOut[Dthread];
  for (int Dcol = 0; Dcol < Nparams; Dcol++) {
    std::vector<double> hessianOutHereCol = hessianOutHere[Dcol];
    for (int Drow = 0; Drow < Nparams; ++Drow) {
      result(Drow, Dcol) += hessianOutHereCol[Drow];
    }
  }
}
return(result);
}



std::vector<double> grad(
  const std::vector<double> parameters,
  const Data& data,
  const Config& config,
  const Rcpp::List& strata,
  const bool onlyRandom
  ) {
  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

  const size_t Nparams = parameters.size();
  const size_t Ngroup = strataP.size()-1;

  std::vector<std::vector<double>> gradOut(config.num_threads);


  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

  #pragma omp parallel
  {
    std::vector<double> gradHere(Nparams, 0);

      #pragma omp for
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {
      auto fun=adFunGroup(
        parameters, data, config, strataI, 
        strataP[Dgroup], strataP[Dgroup+1]);

      auto gradThisGroup = fun.Jacobian(parameters);
      for(size_t D=0;D<Nparams;D++) {
        gradHere[D]+= gradThisGroup[D];
      }

    } // group

#pragma omp single 
    {
      // Q likelihood
      const bool useQ = data.Qdiag.size();
      if(useQ) {
        auto fun = adFunQ(parameters, data, config);
        auto gradQ = fun.Jacobian(parameters);
        for(size_t D=0;D<Nparams;D++) {
          gradHere[D]+= gradQ[D];
        }
      } // useQ
    } // single Q
    const int tid=omp_get_thread_num();
    gradOut[tid] = gradHere;
  }// parallel


  const size_t Nout = onlyRandom?data.Ngamma:Nparams;
  const size_t outOffset = onlyRandom?data.Nbeta:0;

  std::vector<double> result(Nout, 0);

  for(size_t Dthread=0;Dthread<config.num_threads;++Dthread) {
    std::vector<double> outHere = gradOut[Dthread];
    for (int Drow = 0; Drow < Nout; ++Drow) {
      result[Drow] += outHere[Drow+outOffset];
    }
  }

  return(result);

}

/* R exported stuff */

//' @export
// [[Rcpp::export]]
double jointLogDens(
  Rcpp::NumericVector x, 
  Rcpp::List data, 
  Rcpp::List config
  ) {

  Data   dat(data);
  Config cfg(config);
 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(x);

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
  const Rcpp::List& strata,
  const bool onlyRandom
  ) {

 const Data dataC(data);
 const Config configC(config);

 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);

 std::vector<double> result1 = grad(parametersC, dataC, configC, strata, onlyRandom);

 Rcpp::NumericVector result(result1.begin(), result1.end());
 
 return(result);
}


//' @export
// [[Rcpp::export]]
Rcpp::S4 hessian(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config,
  const Rcpp::List& strata,
  const Rcpp::List& sparsity,
  const bool onlyRandom
  ) {

 const Data dataC(data);
 const Config configC(config);

 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);

 Rcpp::S4 result = hessian(
  parametersC,dataC, configC, strata, sparsity, onlyRandom
  );
 return(result);
}


//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hessianQdense(
      const Rcpp::NumericVector parameters, 
      const Rcpp::List& data,
      const Config &config) {

      const Data dataC(data);
      const Config configC(config);

      std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);
      Rcpp::NumericMatrix result = hessianQdense(parametersC, dataC, configC);
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hessianDense(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config,
  const Rcpp::List& strata
  ) {

 const Data dataC(data);
 const Config configC(config);

 std::vector<double> parametersC = Rcpp::as<std::vector<double>>(parameters);

Rcpp::NumericMatrix result = hessianDense(
  parametersC, dataC, configC, strata
  );
 return(result);
}

