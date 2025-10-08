#include"hpol.hpp"
#include"loglikHelpers.hpp"

// #define DEBUG

#include<omp.h>
#include <set>


template<class Type> CppAD::vector<Type> loglikOneStrata(
const int Dstrata,
const CppAD::vector<Type>& gamma,
const PackedParams<double>& parameters, // gamma is ignored
const Data& data, 
const Config& cfg
){

  CppAD::vector<Type> result(1);

  CppAD::vector<Type> etaHere = compute_eta_for_stratum<Type, double>(
    Dstrata, data, gamma, parameters.beta);

  auto contrib = accumulate_contrib_for_stratum<Type, double>(
    Dstrata, data, etaHere, parameters, cfg
    );

  result[0]= contrib;
  return(result);
}


// Compute scaled gamma values and accumulate log-likelihood contribution
template <class TypeGamma, class TypeTheta>
CppAD::vector<TypeGamma>  loglikQ(
    const CppAD::vector<TypeGamma>& gamma,  
    const PackedParams<TypeTheta>& latent,
    const Data& data
) {

    CppAD::vector<TypeGamma> gammaScaled(data.Ngamma);
    CppAD::vector<TypeGamma> result(1);
    result[0] = TypeGamma(0);

    for (size_t D = 0; D < data.Ngamma; ++D) {
        size_t mapHere = data.map[D];

        gammaScaled[D] = gamma[D] / latent.theta[mapHere];

        result[0] += TypeGamma(latent.logTheta[mapHere])
                      + TypeGamma(0.5* data.Qdiag[D]) * gammaScaled[D] * gammaScaled[D] ;
    }

      // Q offdiag    
    for(size_t D = 0; D < data.Nq; D++) {
        result[0] += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] 
          * TypeGamma(data.QsansDiag.x[D]);
    }

    return(result);
}

template<class Type> CppAD::vector<Type> loglikNStrata(
const int Dstrata,
const int Niter,
const CppAD::vector<Type>& gamma,
const PackedParams<double>& parameters, // gamma is ignored
const Data& data, 
const Config& cfg
){
  const size_t end = data.Nstrata < (Dstrata+Niter) ? data.Nstrata : (Dstrata+Niter);
  CppAD::vector<Type> result(1);
  result[0] = Type(0);

  for(size_t Diter=Dstrata; Diter < end; ++Diter ) {
    result[0] += loglikOneStrata(Dstrata, gamma, parameters, data, cfg)[0];
  }
  return(result);
}


//' @export
// [[Rcpp::export]]
double objectiveFunctionNoDiff(
  Rcpp::NumericVector parameters, 
  Rcpp::List dataList, 
  Rcpp::List configList
  ) {

  Data   data(dataList);
  Config cfg(configList);

  auto latent = unpack_params(parameters, data, cfg);

  std::vector<double> loglik(cfg.num_threads);
  CppAD::vector<double> gammaScaled(data.Ngamma);


  omp_set_num_threads(cfg.num_threads);

      #pragma omp parallel
  { 

    const int tid=omp_get_thread_num();
    loglik[tid] = double(0.0);

    #pragma omp for nowait
    for (size_t Dstrata = 0; Dstrata < data.Nstrata; Dstrata++) {

      loglik[tid] += loglikOneStrata<double>(
        Dstrata,
        latent.gamma,
        latent,
        data, 
        cfg
        )[0];
      }

  // Q diag.  
        #pragma omp for nowait
      for(size_t D=0;D<data.Ngamma;D++) {
        size_t mapHere = data.map[D];

        gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
        loglik[tid] += latent.logTheta[mapHere] +
          0.5*gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
      }
  } // end parallel block


  double result=0.0;
  for(size_t D=0;D<loglik.size();D++) {
    result += loglik[D];
  }

  // Q offdiag    
      for(size_t D = 0; D < data.Nq; D++) {
        result += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] * data.QsansDiag.x[D];
    }

  return result;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector objectiveFunctionGrad(
  Rcpp::NumericVector parameters, 
  Rcpp::List dataList, 
  Rcpp::List configList
  ) {

// options: chunks of strata, or build tape once outside loop

  const Data   data(dataList);
  const Config cfg(configList);
  const PackedParams<double> parameters_extra = unpack_params(parameters, data, cfg);

  CppAD::vector<CppAD::vector<double>>  grad_local(cfg.num_threads);
  for(size_t D=0;D<cfg.num_threads;D++) {
    grad_local[D].resize(parameters_extra.Ngamma);
  }

  if(cfg.verbose) Rcpp::Rcout << "Ngamma " << parameters_extra.Ngamma << " starting parallel " << cfg.num_threads << " threads\n";

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
    CppAD::vector<double> w(1); w[0] = 1.0;   

    for (size_t D = 0; D < parameters_extra.Ngamma; D++ ) {
      const double g = parameters_extra.gamma[D];
      ad_paramsT[D] = g;
      x_valT[D] = g;
      grad_localT[D]=0;
    }

   #pragma omp for 
    for (size_t Dstrata = 0; Dstrata < data.Nstrata;
      Dstrata+= cfg.strataPerIter) {


    CppAD::Independent(ad_paramsT);  // Tell CppAD these are inputs for differentiation

    auto yvalT = loglikNStrata<CppAD::AD<double>>(
      Dstrata,
      cfg.strataPerIter,
      ad_paramsT,
      parameters_extra,
      data, 
      cfg
      );

    auto funT = CppAD::ADFun<double>(ad_paramsT, yvalT );
    funT.optimize();
    funT.Forward(0, x_valT);   
    auto grad_resultT = funT.Reverse(1, w); 

    for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
      grad_localT[Dparam] += grad_resultT[Dparam];
    }

  } // Dstrata

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
    grad_total[Dparam] = acc;
  }

// Q 
  CppAD::vector<CppAD::AD<double>> ad_params(parameters_extra.Ngamma);
  CppAD::vector<double> x_val(parameters_extra.Ngamma);

  for (size_t D = 0; D < parameters_extra.Ngamma; D++) {
      const double g = parameters_extra.gamma[D];
      ad_params[D] = g;
      x_val[D] = g;
  }


  CppAD::Independent(ad_params);  
  auto yvec = loglikQ(ad_params,  parameters_extra, data); 

  auto funQ = CppAD::ADFun<double>(ad_params, yvec);
  funQ.optimize();
  auto grad = funQ.Jacobian(x_val);

  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
    grad_total[Dparam] += grad[Dparam];
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

//' @export
// [[Rcpp::export]]
Rcpp::List objectiveFunctionHessian2(
  Rcpp::NumericVector parameters, 
  Rcpp::List dataList, 
  Rcpp::List configList
  ) {


  const Data   data(dataList);
  const Config cfg(configList);
  const PackedParams<double> parameters_extra = unpack_params(parameters, data, cfg);

  const Rcpp::List sparsity = cfg.sparsity;
  const Rcpp::List pairs = sparsity["random"];
  auto pairsi = get_intvec_copy(pairs, "i");
  auto pairsj = get_intvec_copy(pairs, "j");
  auto pairsp = get_intvec_copy(pairs, "p");
  auto jForDiag = get_intvec_copy(pairs, "jNoOffDiag");
  auto Npairs = pairsj.size();
  auto NoffDiag = pairsi.size();


  CppAD::vector<CppAD::vector<double>>  local(cfg.num_threads);
  CppAD::vector<CppAD::vector<double>>  localOffdiag(cfg.num_threads);

  for(size_t D=0;D<cfg.num_threads;D++) {
    local[D].resize(parameters_extra.Ngamma);
    localOffdiag[D].resize(NoffDiag);
  }

  omp_set_num_threads(cfg.num_threads);

  CppAD::thread_alloc::parallel_setup(
    cfg.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );


// diagonals
      #pragma omp parallel
 { 

    const int tid=omp_get_thread_num();

    CppAD::vector<CppAD::AD<double>> ad_paramsT(parameters_extra.Ngamma);
    CppAD::vector<double> x_valT(parameters_extra.Ngamma);
    CppAD::vector<double> localT(parameters_extra.Ngamma);
    CppAD::vector<double> localOffdiagT(NoffDiag);
    CppAD::vector<double> w(1); w[0] = 1.0;   
    CppAD::vector<double> u(parameters_extra.Ngamma, 0.0);

    // pattern for the diagonals
    std::vector<std::set<size_t>> p(parameters_extra.Ngamma);
    std::vector<size_t> row(parameters_extra.Ngamma), col(parameters_extra.Ngamma);

    for(size_t j=0; j<parameters_extra.Ngamma; ++j) {
      p[j].insert(j);
      row[j]=j;
      col[j]=j;
    }
    CppAD::sparse_hessian_work work; 
    CppAD::vector<double> hes(parameters_extra.Ngamma);

    for (size_t D = 0; D < parameters_extra.Ngamma; ++D) {
      const double g = parameters_extra.gamma[D];
      ad_paramsT[D] = g;
      x_valT[D] = g;
      localT[D]=0;
    }
    for (size_t D = 0; D < NoffDiag; ++D) {
      localOffdiagT[D]=0;
    }


   #pragma omp for 
    for (size_t Dstrata = 0; Dstrata < data.Nstrata; Dstrata++) {


    CppAD::Independent(ad_paramsT);  // Tell CppAD these are inputs for differentiation

    auto yvalT = loglikOneStrata<CppAD::AD<double>>(
      Dstrata,
      ad_paramsT,
      parameters_extra,
      data, 
      cfg
      );

    auto funT = CppAD::ADFun<double>(ad_paramsT, yvalT );
    funT.optimize();
    work.clear();
    size_t count=funT.SparseHessian(x_valT, w, p, row, col, hes, work);

    for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
      localT[Dparam] += hes[Dparam];
    }

    // off diagonals
    funT.Forward(0, x_valT);   
    for(size_t Doffdiag=0; Doffdiag < Npairs; ++Doffdiag) {
        u[pairsj[Doffdiag]] = 1;
        funT.Forward(1, u);
        auto resultT = funT.Reverse(2, u); 
        u[pairsj[Doffdiag]] = 0;
        auto pairsEnd = pairsp[Doffdiag+1];
        for(size_t Di=pairsp[Doffdiag];Di<pairsEnd;Di++) {
          localOffdiagT[Di] += resultT[2*pairsi[Di]+1];
        }

    }


  } // Dstrata

   for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
      local[tid][Dparam] = localT[Dparam];
  }
   for(size_t Dparam = 0; Dparam < NoffDiag; Dparam++) {
      localOffdiag[tid][Dparam] = localOffdiagT[Dparam];
  }



} // parallel

 if(cfg.verbose) Rcpp::Rcout << "done parallel\n";


  Rcpp::NumericVector diag_total(parameters_extra.Ngamma);
  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; ++Dparam){
    double acc = 0.0;
    for(size_t Dthread=0; Dthread < cfg.num_threads; ++Dthread) {
      acc += local[Dthread][Dparam];
    }
    diag_total[Dparam] = acc;
  }

  Rcpp::NumericVector off_diag_total(NoffDiag);
  for(size_t Dparam = 0; Dparam < NoffDiag; ++Dparam){
    double acc = 0.0;
    for(size_t Dthread=0; Dthread < cfg.num_threads; ++Dthread) {
      acc += localOffdiag[Dthread][Dparam];
    }
    off_diag_total[Dparam] = acc;
  }

#ifdef UNDEF
// Q 
  CppAD::vector<CppAD::AD<double>> ad_params(parameters_extra.Ngamma);
  CppAD::vector<double> x_val(parameters_extra.Ngamma);
  CppAD::vector<double> w(1); w[0] = 1.0;   

  for (size_t D = 0; D < parameters_extra.Ngamma; D++) {
      const double g = parameters_extra.gamma[D];
      ad_params[D] = g;
      x_val[D] = g;
  }

  CppAD::Independent(ad_params);  
  auto yvec = loglikQ(ad_params,  parameters_extra, data); 

  auto funQ = CppAD::ADFun<double>(ad_params, yvec);
  funQ.optimize();
  funQ.Forward(0, x_val);   
  auto resultQ = funQ.Reverse(2, w); 

  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
    diag_total[Dparam] += resultQ[Dparam];
  }
#endif


  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("diag") = diag_total,
    Rcpp::Named("off_diag") = off_diag_total
    );
  return(result);

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

    loglik += contrib;

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

