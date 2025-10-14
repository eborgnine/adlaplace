
#include"hpol.hpp"
#include"loglikHelpers.hpp"




// #define DEBUG

#include<omp.h>



//' @export
// [[Rcpp::export]]
double objectiveFunctionNoDiff(
  Rcpp::NumericVector x, 
  Rcpp::List data, 
  Rcpp::List config
  ) {

  Data   dat(data);
  Config cfg(config);

  auto latent = unpack_params(x, dat, cfg);

  std::vector<double> loglik(cfg.num_threads);
  std::vector<double> Qpart(cfg.num_threads);
  CppAD::vector<double> gammaScaled(dat.Ngamma);


  omp_set_num_threads(cfg.num_threads);

      #pragma omp parallel
  { 

    const int tid=omp_get_thread_num();
    double loglikT = double(0.0), QpartT=double(0.0);

    #pragma omp for nowait
    for (size_t Dstrata = 0; Dstrata < dat.Nstrata; Dstrata++) {
      loglikT += loglikOneStrata<double>(
        Dstrata,
        latent.gamma,
        latent,
        dat, 
        cfg
        )[0];
      }

  // Q diag.  
        #pragma omp for nowait
      for(size_t D=0;D<dat.Ngamma;D++) {
        size_t mapHere = dat.map[D];

        gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
        QpartT += latent.logTheta[mapHere] +
        0.5*gammaScaled[D]*gammaScaled[D]*dat.Qdiag[D];
      }
      loglik[tid] = loglikT;      
      Qpart[tid] = QpartT;
  } // end parallel block


  double resultL=0.0, resultQ = 0.0;
  for(size_t D=0;D<loglik.size();D++) {
    resultL += loglik[D];
    resultQ += Qpart[D];
  }
  double result = -resultL + resultQ;




  // Q offdiag    
  for(size_t D = 0; D < dat.Nq; D++) {
    result += gammaScaled[dat.QsansDiag.i[D]] * gammaScaled[dat.QsansDiag.j[D]] * dat.QsansDiag.x[D];
  }
  if (cfg.verbose ) Rcpp::Rcout << "L " << resultL << " Q " << resultQ << 
    " total " << result << "\n";

  return result;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector objectiveFunctionGrad(
  Rcpp::NumericVector x, 
  Rcpp::List data, 
  Rcpp::List config
  ) {

// options: chunks of strata, or build tape once outside loop

  const Data   dat(data);
  const Config cfg(config);
  const PackedParams<double> parameters_extra = unpack_params(x, dat, cfg);
  const size_t strataPerIter = cfg.strataPerIter > 0? cfg.strataPerIter : 1;

  CppAD::vector<CppAD::vector<double>>  grad_local(cfg.num_threads);

  for(size_t D=0;D<cfg.num_threads;D++) {
    grad_local[D].resize(parameters_extra.Ngamma);
  }

  if(cfg.verbose) Rcpp::Rcout << "Ngamma " << parameters_extra.Ngamma << " Nstrata " << 
      dat.Nstrata << " starting parallel " << cfg.num_threads << " threads\n";

  omp_set_num_threads(cfg.num_threads);

  CppAD::thread_alloc::parallel_setup(
    cfg.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );

      #pragma omp parallel
  {

    const int tid=omp_get_thread_num();

    CppAD::vector<CppAD::AD<double>> ad_paramsT(parameters_extra.Ngamma);
    CppAD::vector<double> x_valT(parameters_extra.Ngamma);
    CppAD::vector<double> grad_localT(parameters_extra.Ngamma);
    CppAD::vector<double> wT({1.0});   


    for (size_t D = 0; D < parameters_extra.Ngamma; D++ ) {
      const double g = parameters_extra.gamma[D];
      ad_paramsT[D] = g;
      x_valT[D] = g;
      grad_localT[D]=0;
    }

   #pragma omp for nowait
    for (size_t DstrataOuter = 0; 
      DstrataOuter < dat.Nstrata; 
      DstrataOuter += cfg.strataPerIter
      )     {

    CppAD::Independent(ad_paramsT);  // Tell CppAD these are inputs for differentiation

    auto yvalT = loglikNStrata<CppAD::AD<double>>(
      DstrataOuter,
      cfg.strataPerIter,
      ad_paramsT,
      parameters_extra,
      dat, 
      cfg
      );

    CppAD::ADFun<double> funT(ad_paramsT, yvalT);
    CppAD::vector<double> y0 = funT.Forward(0, x_valT);
    auto grad_resultT = funT.Reverse(1, wT); 

      for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
        grad_localT[Dparam] += grad_resultT[Dparam];
      }
  } // DstrataOuter



  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
    grad_local[tid][Dparam] = grad_localT[Dparam];
  }

  } // end parallel block

  if(cfg.verbose) Rcpp::Rcout << "done parallel\n";

  Rcpp::NumericVector grad_total(parameters_extra.Ngamma);
  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; ++Dparam){
    double acc = 0.0;
    for(size_t Dthread=0; Dthread < cfg.num_threads; ++Dthread) {
      acc += grad_local[Dthread][Dparam];
    }
    grad_total[Dparam] = - acc;
  }


  if(cfg.verbose) Rcpp::Rcout << "gl " << grad_local[0][0] <<  " d0 " << grad_total[0] << " d1 " << grad_total[1] << "\n";

// Q , put this in parallel?
  CppAD::vector<CppAD::AD<double>> ad_params(parameters_extra.Ngamma);
  CppAD::vector<double> x_val(parameters_extra.Ngamma);
  CppAD::vector<double> w({1.0});   

  for (size_t D = 0; D < parameters_extra.Ngamma; D++) {
    const double g = parameters_extra.gamma[D];
    ad_params[D] = g;
    x_val[D] = g;
  }


  CppAD::Independent(ad_params);  
  auto yvec = loglikQ(ad_params,  parameters_extra, dat);
  auto funQ = CppAD::ADFun<double>(ad_params, yvec);
  funQ.Forward(0, x_val);   
  auto gradQ= funQ.Reverse(1, w); 

  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
    grad_total[Dparam] += gradQ[Dparam];
  }

  return grad_total;
}



//' @export
// [[Rcpp::export]]
Rcpp::List objectiveFunctionHessian(
  Rcpp::NumericVector parameters, 
  Rcpp::List dataList, 
  Rcpp::List configList
  ) {

  Data   data(dataList);
  Config config(configList);
  const PackedParams<double> parameters_extra = unpack_params(parameters, data, config);

  const Rcpp::List sparsity = config.sparsity;
  const Rcpp::List pairs = sparsity["random"];
  auto pairsi = get_intvec_copy(pairs, "i");
  auto pairsj = get_intvec_copy(pairs, "j");
  auto pairsp = get_intvec_copy(pairs, "p");
  auto onesJ = get_intvec_copy(pairs, "onesJ");
  auto onesI = get_intvec_copy(pairs, "onesI");
  auto Npairs = pairsj.size();
  auto NoffDiag = pairsi.size();
  auto Nones = onesJ.size();

  auto Nparams = parameters_extra.Ngamma;
  auto x_val = parameters_extra.gamma;

  std::vector<double> resultOffDiag(NoffDiag * ( config.dense ? Nparams : 1 ) );
  std::vector<double> resultDiag(Nones);
  std::vector<double> resultColSum, diagTest;

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );


  if (config.verbose ) Rcpp::Rcout << "starting parallel " << resultOffDiag.size() << " size " << 
    Nones << " diags " << config.num_threads << " threads\n";



  #pragma omp parallel
  {


    std::vector<double> y_val(1);
    std::vector<double> w{0, 1.0};  
    std::vector<double> direction(Nparams, 0.0);
    std::vector<double> directionZeros(Nparams, 0.0);
    std::vector<double> forward2direction(3*Nparams, 0.0);


    CppAD::vector<CppAD::AD<double>> ad_params(Nparams);  
    for (size_t D = 0; D < Nparams; D++) {
      ad_params[D] = x_val[D];  // Initialize CppAD variables
      forward2direction[3*D] = x_val[D];
    }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation

  auto y = objectiveFunctionInternal(ad_params, data, config);  
  CppAD::ADFun<double> fun(ad_params, y);
  fun.Forward(0, x_val);


#pragma omp single
  {
    std::fill(direction.begin(), direction.end(), 1.0);
    fun.Forward(0, x_val);
    fun.Forward(1, direction);
    resultColSum = fun.Reverse(2, w);

    size_t Dj = onesJ[0];
    forward2direction[3*Dj + 1] = 1.0;
    diagTest = fun.Forward(2, forward2direction); 
    forward2direction[3*Dj + 1] = 0.0;
  }


    #pragma omp for nowait
  for(int Dpair=0; Dpair < Npairs; ++Dpair) {
    auto Dj = pairsj[Dpair];

    std::fill(direction.begin(), direction.end(), 0.0);
    direction[Dj]  = 1.0;     

    fun.Forward(1, direction);
    auto resRev = fun.Reverse(2, w);

    if(config.dense) {
      for(size_t Di=0,Dindex=(Dpair*Nparams);Di<Nparams;Di++,Dindex++) {
        resultOffDiag[Dindex] = resRev[2*Di];
      }
    } else {
      auto IindexEnd = pairsp[Dpair+1];
      for(auto Diindex=pairsp[Dpair]; Diindex < IindexEnd; Diindex++) {
        auto Di = pairsi[Diindex];
        resultOffDiag[Diindex] = resRev[2*Di];
      }
    }

  } // Dpair


  // diagonals with forward.  do this for columns with only one off-diagonal?
    #pragma omp for nowait
  for(size_t Ddiag=0; Ddiag < Nones; ++Ddiag) {
    size_t Dj = onesJ[Ddiag];
    forward2direction[3*Dj + 1] = 1.0;
    auto   resForward = fun.Forward(2, forward2direction); 
    resultDiag[Ddiag] = 2*resForward[2];
    forward2direction[3*Dj + 1] = 0.0;
  }

} //parallel

if (config.verbose ) Rcpp::Rcout << "done parallel\n";


Rcpp::List resultList = Rcpp::List::create(
  Rcpp::Named("diag") = Rcpp::wrap(resultDiag),
  Rcpp::Named("off_diag") = Rcpp::wrap(resultOffDiag),
  Rcpp::Named("colSum") = Rcpp::wrap(resultColSum),
  Rcpp::Named("test") = Rcpp::wrap(diagTest)
  );
return(resultList);
}



template<class Type>
CppAD::vector<Type>  objectiveFunctionInternal(
 const CppAD::vector<Type>& ad_params,  
 const Data& data,
 const Config& config
 ) {

  auto latent=unpack_params<Type>(ad_params, data, config);

  CppAD::vector<Type> minusLogDens(1);
  Type loglik = Type(0);


  for (size_t Dstrata = 0; Dstrata < data.Nstrata;  Dstrata++) {

    auto etaHere = compute_eta_for_stratum(
      Dstrata, data, latent.gamma, latent.beta);

    auto contrib = accumulate_contrib_for_stratum<Type, Type>(
      Dstrata, data, etaHere, latent, config);

    loglik += contrib[0];

  }

  auto randomContribution = loglikQ(latent.gamma, latent, data);



  minusLogDens[0] =  - loglik + randomContribution[0];


#ifdef EVALCONSTANTS
  minusLogDens[0] += Ngamma * HALFLOGTWOPI;
  minusLogDens[0] -= Rcpp::as<double>(config.halfLogDetQ);
#endif


  return minusLogDens;
}



//' @export
// [[Rcpp::export]]
Rcpp::List objectiveFunctionC(
  Rcpp::NumericVector parameters, 
  Rcpp::List dataList, 
  Rcpp::List configList
  ) {

  Data   data(dataList);
  Config config(configList);

  size_t Nparams = parameters.size();


  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  
  
  auto y = objectiveFunctionInternal<CppAD::AD<double>>(ad_params, data, config);

  if(config.verbose ) {
    Rcpp::Rcout << "y " << y[0] << "\n";
  }
  CppAD::ADFun<double> fun(ad_params, y);

  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) { 
    x_val[i] = parameters[i];
  }

  if (config.verbose ) {
    Rcpp::Rcout << "forward: ";
  }

  std::vector<double> y_val(1);
  y_val = fun.Forward(0, x_val);

  if (config.verbose ) {
    Rcpp::Rcout << "f0 " << y_val[0] << "\n";
  }


  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("value") = y_val[0]
    );

  if(config.maxDeriv == 0) {
    return result;
  }

  if (config.verbose ) {
    Rcpp::Rcout << "grad ";
  }
  std::vector<double> grad = fun.Jacobian(x_val);
  if (config.verbose ) {
    Rcpp::Rcout << ".\n";
  }
  result["grad"] = grad;
  if(config.maxDeriv == 1) {
    return result;
  }

// hessian
  if (config.verbose ) {
    Rcpp::Rcout << "hess " << config.num_threads << " threads\n";
  }



  bool dense = config.dense;
  if(!config.sparsity.size()) {
    Rcpp::warning("no sparsity pattern specified, producing dense matrix");
    dense = true;
  }

  Rcpp::NumericMatrix denseHessianOut;
  Rcpp::IntegerVector Hrow, Hcol, Hp;
  Rcpp::NumericVector Hvalue;
  Rcpp::List sparsityR, sparsitySecond;

  if(dense) {
    denseHessianOut = Rcpp::NumericMatrix(Nparams, Nparams);
  } else {
    sparsitySecond = config.sparsity["second"];
    if(config.theta.size()) {
      // hessian only for gammas
      if (config.verbose ) {
        Rcpp::Rcout << "sparsity only for random effects\n";
      }   
      sparsityR = sparsitySecond["random"];      
    } else {
      if (config.verbose ) {
        Rcpp::Rcout << "full parameter sparsity\n";
      }   
      sparsityR = sparsitySecond["full"];      
    }
    Hrow = sparsityR["i"]; 
    Hcol = sparsityR["j"];
    Hp = sparsityR["p"];
    Hvalue = Rcpp::NumericVector(Hrow.size());
  }

  omp_set_num_threads(config.num_threads);

    // Replicate fun object for each thread
  std::vector<CppAD::ADFun<double>> fun_threads(config.num_threads);
  for (int i = 0; i < config.num_threads; ++i) {
    fun_threads[i] = fun;
  }

  const std::vector<double> w(1, 1.0);

    #pragma omp parallel
  {

    const int tid = omp_get_thread_num();
    std::vector<double> u(Nparams, 0.0);

    #pragma omp for
    for (int Dcol = 0; Dcol < Nparams; Dcol++) {
      std::fill(u.begin(), u.end(), 0.0);
      u[Dcol] = 1.0;
      fun_threads[tid].Forward(0, x_val);
      fun_threads[tid].Forward(1, u);
      std::vector<double> ddw = fun_threads[tid].Reverse(2, w);

      if(config.dense) {
        for (int Drow = 0; Drow < Nparams; ++Drow) {
          denseHessianOut(Drow, Dcol) = ddw[2 * Drow + 1];
        }
      } else {
        int HpDcolp1 = Hp[Dcol+1];
        for (int Dindex = Hp[Dcol]; Dindex < HpDcolp1; ++Dindex) {
          int Drow = Hrow[Dindex];
          Hvalue[Dindex] = ddw[2 * Drow + 1];
        }
          } // sprse
      } //  column


  } // parallel


  Rcpp::S4 hessianR = make_CMatrix(
    Hvalue, Hrow, Hp); 

  result["hessian"] = hessianR;

  if(config.dense) {
    result["denseHessian"] = denseHessianOut;
  }

  if (config.verbose ) { 
    Rcpp::Rcout << " done." << std::endl;
  }

  return result;
}

