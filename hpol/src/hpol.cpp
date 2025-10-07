#include"hpol.hpp"
#include"loglikHelpers.hpp"

#define DEBUG

#include<omp.h>



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

  result[0]= Type(1);//etaHere[0]; //contrib;

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

  omp_set_num_threads(cfg.num_threads);

  CppAD::vector<double> loglik(cfg.num_threads);
  CppAD::vector<int> NperThread(cfg.num_threads);
  CppAD::vector<double> gammaScaled(data.Ngamma);

  CppAD::thread_alloc::hold_memory(true);

if(cfg.verbose) {
    Rcpp::Rcout << "t " << cfg.num_threads << " s " << data.Nstrata << " gamma0 " << latent.gamma[0] << "\n";
}


      #pragma omp parallel
  { 

    const int tid=omp_get_thread_num();
    loglik[tid] = double(0.0);
    NperThread[tid] = 0;


    #pragma omp for nowait
    for (size_t Dstrata = 0; Dstrata < 800;//data.Nstrata; 
      Dstrata++) {

      loglik[tid] += loglikOneStrata<double>(
        Dstrata,
        latent.gamma,
        latent,
        data, 
        cfg
        )[0];
      NperThread[tid]++; 
      }
#ifdef UNDEF
  // Q diag.  
        #pragma omp for nowait
      for(size_t D=0;D<data.Ngamma;D++) {
        size_t mapHere = data.map[D];

        gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
        loglik[tid] += latent.logTheta[mapHere] +
          0.5*gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
      }
#endif
  } // end parallel block

if(cfg.verbose) {
  for(size_t D=0;D<loglik.size();D++) {
    Rcpp::Rcout << "t " << D << " n " << NperThread[D] << " l " << loglik[D] << "\n";
  }
}

  double result=0.0;
  for(size_t D=0;D<loglik.size();D++) {
    result += loglik[D];
  }
#ifdef UNDEF
  // Q offdiag    
      for(size_t D = 0; D < data.Nq; D++) {
        result += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] * data.QsansDiag.x[D];
    }
#endif

  return result;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector objectiveFunctionGrad(
  Rcpp::NumericVector parameters, 
  Rcpp::List dataList, 
  Rcpp::List configList
  ) {

  const Data   data(dataList);
  const Config cfg(configList);
  const PackedParams<double> parameters_extra = unpack_params(parameters, data, cfg);

  omp_set_num_threads(cfg.num_threads);

  CppAD::vector<CppAD::vector<CppAD::AD<double>>> ad_params(cfg.num_threads);
  CppAD::vector<CppAD::vector<double>>  grad_local(cfg.num_threads);
  CppAD::vector<CppAD::ADFun<double>> fun(cfg.num_threads);
  std::vector<std::vector<double> > x_val(cfg.num_threads);

  if(cfg.verbose) Rcpp::Rcout << "Ngamma " << parameters_extra.Ngamma << " starting parallel " << cfg.num_threads << " threads\n";


  CppAD::thread_alloc::hold_memory(true);
      #pragma omp parallel
  { 

    const int tid=omp_get_thread_num();


    grad_local[tid] = CppAD::vector<double>(parameters_extra.Ngamma);
    ad_params[tid] = CppAD::vector<CppAD::AD<double>>(parameters_extra.Ngamma);
    x_val[tid] = std::vector<double>(parameters_extra.Ngamma);

    for (size_t D = 0; D < parameters_extra.Ngamma; D++) {
      grad_local[tid][D] =0.0;
      ad_params[tid][D] = x_val[tid][D] = parameters_extra.gamma[D];  
    }

    CppAD::Independent(ad_params[tid]);  // Tell CppAD these are inputs for differentiation


//   #pragma omp for nowait
//    for (size_t Dstrata = 0; Dstrata < data.Nstrata; Dstrata++) {
      size_t Dstrata = tid;
      auto yvec = loglikOneStrata<CppAD::AD<double>>(
        Dstrata,
        ad_params[tid],
        parameters_extra,
        data, 
        cfg
        );

      fun[tid] = CppAD::ADFun<double>(ad_params[tid], yvec);
      auto grad = fun[tid].Jacobian(x_val[tid]);

      for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
        grad_local[tid][Dparam] += grad[Dparam];
      }

//  } // Dstrata


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

#ifdef UNDEF
// Q 

 for (size_t D = 0; D < parameters_extra.Ngamma; D++) {
      ad_params[0][D] = x_val[0][D] = parameters_extra.gamma[D];  
  }

  auto yvec = loglikQ(ad_params[0],  parameters_extra, data); 

  auto funQ = CppAD::ADFun<double>(ad_params[0], yvec);
  auto grad = funQ.Jacobian(x_val[0]);

  for(size_t Dparam = 0; Dparam < parameters_extra.Ngamma; Dparam++) {
    grad_total[Dparam] += grad[Dparam];
  }
#endif

  return grad_total;
}



template<class Type>
CppAD::vector<Type>  objectiveFunctionInternal(
 const CppAD::vector<Type>& ad_params,  
 const Data& data,
 const Config& config
 ) {

  Rcpp::Rcout << "a\n";
  auto latent=unpack_params<Type>(ad_params, data, config);
  Rcpp::Rcout << "b\n";

  CppAD::vector<Type> minusLogDens(1);
  minusLogDens[0]=0.0;

#ifdef DEBUG
  Rcpp::Rcout << " logSqrtNU " << latent.logSqrtNu << 
  " oneOverSqrtNu " << latent.oneOverSqrtNu << " lgammaOneOverSqrtNu " <<
  latent.lgammaOneOverSqrtNu << std::endl << "theta ";
  for(int Dtheta =0;Dtheta < latent.theta.size();++Dtheta) {
    Rcpp::Rcout << Dtheta << " " << latent.theta[Dtheta] << " " << latent.logTheta[Dtheta] << std::endl;
  }
#endif    

  omp_set_num_threads(config.num_threads);

  CppAD::vector<Type>  loglik(config.num_threads);


      #pragma omp parallel
  { 

    // elements beta, gamma, theta, logtheta


    const int tid=omp_get_thread_num();
    loglik[tid] = Type(0);


    #pragma omp for nowait
    for (size_t Dstrata = 0; Dstrata < data.Nstrata; 
      Dstrata++) {

      CppAD::vector<Type> etaHere = compute_eta_for_stratum<Type, Type>(
        Dstrata, data, latent.gamma, latent.beta);

    auto contrib = accumulate_contrib_for_stratum<Type, Type>(
      Dstrata, data, etaHere, latent, config
      );

    loglik[tid] += contrib;

  }

} //parellel block


  // random effect density
    // log(|Q|) + 0.5 * gamma^T Q gamma
  // diagonals
CppAD::vector<Type> gammaScaled(data.Ngamma);

Type randomContributionDiag = Type(0);
Type offdiagQ               = Type(0);

for(size_t D=0;D<data.Ngamma;D++) {
  size_t mapHere = data.map[D];

  gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
  randomContributionDiag += latent.logTheta[mapHere] +
  0.5*gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
}

    // Q offdiag    
for(size_t D = 0; D < data.Nq; D++) {
  offdiagQ += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] * data.QsansDiag.x[D];
}


Type  loglikS = Type(0);



for (int t = 0; t < config.num_threads; ++t) {
#ifdef DEBUG
  Rcpp::Rcout << "thread " << t << "logLik " << loglik[t] << " offdiag " << offdiagQ <<
  " diag " << randomContributionDiag << std::endl;
#endif 
  loglikS += loglik[t];
}

#ifdef DEBUG
Rcpp::Rcout << "logLik " << loglikS << " offdiag " << offdiagQ <<
" diag " << randomContributionDiag << std::endl;
#endif 


minusLogDens[0] =  - loglikS + offdiagQ  + randomContributionDiag;
//  etaLogSum[0]  + etaLogSum[1]


//    Rcpp::Rcout << " " << etaLogSum[0] << " " << etaLogSum[1] << "\n";

#ifdef EVALCONSTANTS
minusLogDens[0] += Ngamma * HALFLOGTWOPI;
minusLogDens[0] -= Rcpp::as<double>(config.halfLogDetQ);
#endif

#ifdef DEBUG
Rcpp::Rcout << "all done " << minusLogDens[0] << std::endl;
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

