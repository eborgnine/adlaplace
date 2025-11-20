
#include"hpol.hpp"

#define DEBUG


/*#if defined(CPPAD_HAS_COLPACK) && CPPAD_HAS_COLPACK
static const char* JAC_COLOR = "colpack";
static const char* HESS_COLOR = "colpack.symmetric";
#else*/
  static const char* JAC_COLOR = "cppad";  // fallback if ColPack not available
  static const char* HESS_COLOR = "cppad.symmetric";
//#endif




// to do: put Q in adpack

  double jointLogDens(
    CppAD::vector<double> parameters, 
    Data& data, 
    Config& config
    ) {
    auto latent = unpack_params(parameters, data, config);

    std::vector<double> logLik(config.num_threads);
    std::vector<double> Qpart(config.num_threads);
    CppAD::vector<double> gammaScaled(data.Ngamma);

    const size_t NqDiag = data.Qdiag.size();

    double logLikSum = 0.0;
    double QpartSum  = 0.0;

    omp_set_num_threads(config.num_threads);

      #pragma omp parallel reduction(+:logLikSum, QpartSum)
    { 

    #pragma omp for nowait
      for (size_t Dstrata = 0; Dstrata < data.Nstrata; Dstrata++) {

        auto etaHere = compute_eta_for_stratum<double, double>(
          Dstrata, data, latent.gamma, latent.beta);
        auto contrib = accumulate_contrib_for_stratum<double, double>(
          Dstrata, data, etaHere, latent, config);

#ifdef DEBUG        
        if(config.verbose) {
          if(CppAD::isnan(contrib[0])) {        
            Rcpp::Rcout << "_" << Dstrata << "_";
          }
        }
#endif

        logLikSum += contrib[0];

      }

  // Q diag.   map could be smaller than gamma.  the first Ngamma-Nmap entries don't have a theta
      if(NqDiag) {
        int NgammaNoMap = data.Ngamma - data.Nmap;
      # pragma omp for
        for(int Dgamma = 0;Dgamma<NgammaNoMap;Dgamma++) {
          gammaScaled[Dgamma] = latent.gamma[Dgamma];
          QpartSum += gammaScaled[Dgamma]*gammaScaled[Dgamma]*(0.5*data.Qdiag[Dgamma]);
        }
      #pragma omp for
        for(int Dmap = 0;Dmap<data.Nmap;Dmap++) {
          int Dgamma=NgammaNoMap + Dmap;
          size_t mapHere = data.map[Dmap];
          gammaScaled[Dgamma] = latent.gamma[Dgamma] / latent.theta[mapHere];
          QpartSum += latent.logTheta[mapHere] + gammaScaled[Dgamma]*gammaScaled[Dgamma]*(0.5*data.Qdiag[Dgamma]);
        }
    }// if NqDiag
  } // end parallel block


  // Q offdiag  Q = Qdiag + Qnd, gamma^T (Qdiag + Qnd) gamma
  double Qoff = 0.0;
  #pragma omp parallel for reduction(+:Qoff) num_threads(config.num_threads)
  for(size_t D = 0; D < data.Nq; D++) {
    Qoff += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] * data.QsansDiag.x[D];
  }

  const double result = -logLikSum + QpartSum + Qoff;

  if (config.verbose ) {
    Rcpp::Rcout << "L " << logLikSum << " Q " << QpartSum << " Qoff " << Qoff <<
    " total " << result << "\n";
  }

  return result;
}





Rcpp::NumericMatrix hessianQdense(
  const CppAD::vector<double> parameters, 
  const Data& data,
  const Config& config) {


  const size_t Nparams = parameters.size();
  Rcpp::NumericMatrix hessianOut(Nparams, Nparams);
  auto fun = adFunQ(parameters, data, config);
  fun.Forward(0, parameters);

  CppAD::vector<double> w(1);
  w[0] = 1.0;
  CppAD::vector<double> u(Nparams);
  for (size_t i = 0; i < Nparams; i++)
    u[i] = 0.0;

  for (int Dcol = 0; Dcol < Nparams; Dcol++) {
    u[Dcol] = 1.0;
    fun.Forward(1, u);
    u[Dcol] = 0.0;

    const auto ddw = fun.Reverse(2, w);
    for (int Drow = 0; Drow < Nparams; ++Drow) {
      hessianOut(Drow, Dcol) = ddw[2 * Drow + 1];
    }
  } // Dcol
  return(hessianOut);
}



GroupPack hessianQsparse(
  const CppAD::vector<double>& parameters, 
  const Data& data, 
  const Config& config
  ) {

  GroupPack Qpack = getAdFunQ(parameters, data, config);

  CppAD::vector<double> w(1);
  w[0] = 1.0;

  if(config.verbose){
    Rcpp::Rcout << "q, colour " << HESS_COLOR << "\n";
  }


  Qpack.fun.sparse_hes(parameters, w,
    Qpack.out_hess, Qpack.pattern_hess, 
    HESS_COLOR, Qpack.work_hess);

  return Qpack;
}


Rcpp::RObject hessian(
  const CppAD::vector<double>& parameters,
  std::vector<GroupPack>& adpack, 
  const Data& data,
  const Config& config
  ) {

  const Rcpp::List sparsity = config.group_sparsity;
  const size_t Ngroup = adpack.size();


  const int N = adpack[0].out_hess.nr();
  const Rcpp::List secondAll = config.sparsity["second"];
  const Rcpp::List secondRandom = secondAll["random"];
  const Rcpp::List secondFull = secondAll["full"];
  const Rcpp::IntegerVector pRandom = secondRandom["p"];
  const Rcpp::IntegerVector pFull= secondFull["p"];
  const int Nrandom = pRandom.size()-1, Nfull = pFull.size()-1;
  const bool onlyRandom = Nrandom == N;

  const Rcpp::IntegerVector iOut = onlyRandom?secondRandom["i"]:secondFull["i"];
  const Rcpp::IntegerVector jOut = onlyRandom?secondRandom["j"]:secondFull["j"];
  const int Nnonzero = iOut.size();
  std::vector<double> hessianOut(Nnonzero);

  if(config.verbose){
    Rcpp::Rcout << "hessian colour " << HESS_COLOR << "\n";
  }

  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
    );


  #pragma omp parallel
  {

    CppAD::vector<double> w(1);
    w[0] = 1.0;
    std::vector<double> outHereAll(Nnonzero, 0.0);

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
      for(size_t D=0;D<Nnonzero;D++) {
        hessianOut[D] +=  outHereAll[D];
      }  
    } // critical

  } //parallel
    // to do: make Q ad fun part of adpack
  auto qRes = hessianQsparse(parameters, data, config);
  const std::vector<size_t>& matchHere = qRes.outRowCol[2];
  const size_t Nhere = matchHere.size();
  const CppAD::vector<double>& hessianOutHere = qRes.out_hess.val();

    for(size_t D=0;D < Nhere; D++) {
      const size_t indexHere = matchHere[D];
      hessianOut[indexHere] += hessianOutHere[D];
    }

    auto result = make_convert_gCmatrix(hessianOut, iOut, jOut, N);


    return(result);
  }



  Rcpp::RObject hessian(
    const CppAD::vector<double>& parameters,
    const Data& data,
    const Config& config
    ) {

    const Rcpp::List strata = config.groups;
    const Rcpp::List sparsity = config.group_sparsity;

    std::vector<GroupPack> adpack = getAdFun(
      parameters, data, config);

    Rcpp::RObject result = hessian(
      parameters, adpack, data, config);
    return(result);
  }

  Rcpp::LogicalMatrix hessianDenseLogical(
    const CppAD::vector<double> parameters,
    const Data& data,
    const Config& config
    ) {
  // single threaded, no Q
    using CppAD::AD;

    const Rcpp::List strata = config.groups;
    const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

    const size_t Nparams = parameters.size();
    const size_t Ngroup  = static_cast<size_t>(strataP.size() - 1);

  // full logical Hessian (we'll OR contributions from each group)
    Rcpp::LogicalMatrix result(Nparams, Nparams);

    if (config.verbose) {
      Rcpp::Rcout << "getting structural Hessian sparsity via for_hes_sparsity\n";
    }

  // loop over groups
    for (size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      if (config.verbose) {
        Rcpp::Rcout << "  group " << Dgroup <<  " N " << Nparams;
      }

    // build ADFun for this group's contribution
      CppAD::ADFun<double> fun = adFunGroup(
        parameters, data, config,
        strataI,
        strataP[Dgroup],
        strataP[Dgroup + 1]
        );

      const size_t n = fun.Domain();
      if (n != Nparams) {
        Rcpp::stop("hessianDenseLogical: fun.Domain() != Nparams");
      }

    // select all domain variables and the single range component
      CPPAD_TESTVECTOR(bool) select_domain(n), select_range(1);
      for (size_t j = 0; j < n; ++j) {
        select_domain[j] = true;
      }
    select_range[0] = true; // scalar output

    // structural Hessian sparsity pattern for this group
    bool internal_bool = true;
    CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t) > pattern;
    fun.for_hes_sparsity(
      select_domain,
      select_range,
      internal_bool,
      pattern
      );

    const CppAD::vector<size_t>& row = pattern.row();
    const CppAD::vector<size_t>& col = pattern.col();
    const size_t K = pattern.nnz();

    // mark all structurally nonzero entries as TRUE;
    // Hessian is symmetric, so set both (i,j) and (j,i)
    if (config.verbose) {
      Rcpp::Rcout << "  size " << K << "\n";
    }

    for (size_t k = 0; k < K; ++k) {
      size_t i = row[k];
      size_t j = col[k];
      result(i, j) = true;
      result(j, i) = true;
    }
  } // Dgroup

  return result;
}


Rcpp::NumericMatrix hessianDense(
  const CppAD::vector<double>& parameters,
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

    CppAD::vector<double> w(1);
    w[0] = 1.0;
    CppAD::vector<double> u(Nparams);
    for (size_t i = 0; i < Nparams; i++)
      u[i] = 0.0;

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
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        u[Dcol] = 0.0;

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

      CppAD::vector<double> w(1);
      w[0] = 1.0;
      CppAD::vector<double> u(Nparams);
      for (size_t i = 0; i < Nparams; i++)
        u[i] = 0.0;

      for (size_t Dcol = 0; Dcol < Nparams; Dcol++) {
        u[Dcol] = 1.0;
        fun.Forward(1, u);
        u[Dcol] = 0.0;
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
  const CppAD::vector<double> parameters,
  std::vector<GroupPack>& adpack, 
  const Data& data,
  const Config& config
  ) {



  const size_t Nparams = parameters.size();
  const size_t Ngroup = adpack.size();
  const bool useQ = data.Qdiag.size()>0;

  if(config.verbose){
    Rcpp::Rcout << "grad, colour " << JAC_COLOR << "\n";
  }

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

    CppAD::vector<double> x = parameters;
    std::vector<double> gradHere(Nparams, 0);
    CppAD::vector<double> w(1);  
    w[0] = 1.0;

      #pragma omp for schedule(dynamic,1) nowait
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

    double* out = gradOut.data();
    const double* here = gradHere.data();
    int n = static_cast<int>(Nparams);

  #  pragma omp critical 
  // pragma parallel reduction(+: out[:n])
    {
      for (int d = 0; d < n; ++d) {
        out[d] += here[d];
      }
    }

  }// parallel

  return(gradOut);
}


std::vector<double> gradDense(
  const CppAD::vector<double> parameters,
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
    CppAD::vector<double> w(1);
    w[0] = 1.0;

      #pragma omp for schedule(dynamic,1) nowait
    for(size_t Dgroup = 0; Dgroup < Ngroup; ++Dgroup) {

      adpack[Dgroup].fun.Forward(0, parameters);
      auto gradThisGroup = adpack[Dgroup].fun.Reverse(1, w);

      for(size_t D=0;D<Nparams;D++) {
        gradHere[D] += gradThisGroup[D];
      }

    } // group


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

    double* out = gradOut.data();
    const double* here = gradHere.data();
    int n = static_cast<int>(Nparams);

    // Every thread runs this whole loop; OpenMP reduces elementwise
  //parallel reduction(+: out[:n])
    #pragma omp critical 
    {
      for (int d = 0; d < n; ++d) out[d] += here[d];
    }

  }// parallel

  return(gradOut);
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
  const size_t N = x.size();
  CppAD::vector<double> parametersC(N);
  for(size_t D=0; D<N;D++) {
    parametersC[D] = x[D];
  }


  std::vector<GroupPack> adpack = getAdFun(
    parametersC, dataC, configC);
  if(configC.verbose) {
    Rcpp::Rcout << "adfun size " << adpack.size() << "\n";
  }

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
  const size_t N = parameters.size();
  CppAD::vector<double> parametersC(N);
  for(size_t D=0; D<N;D++) {
    parametersC[D] = parameters[D];
  }



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

  const size_t N = parameters.size();
  CppAD::vector<double> parametersC(N);
  for(size_t D=0; D<N;D++) {
    parametersC[D] = parameters[D];
  }
  const Data   dataC(data);
  const Config configC(config);

  AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
  std::vector<GroupPack>* packsPtr = ad.ptr;

  // if configC.dense
  std::vector<double> result =configC.dense?
  gradDense(parametersC, *packsPtr, dataC, configC)
  :
  grad(parametersC, *packsPtr, dataC, configC);

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

AdpackHandle ad = getAdpackFromR(adFun, parametersC, dataC, configC);
std::vector<GroupPack>* packsPtr = ad.ptr;


Rcpp::RObject result = hessian(
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

  const size_t N = parameters.size();
  CppAD::vector<double> parametersC(N);
  for(size_t D=0; D<N;D++) {
    parametersC[D] = parameters[D];
  }

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

 const size_t N = parameters.size();
 CppAD::vector<double> parametersC(N);
 for(size_t D=0; D<N;D++) {
  parametersC[D] = parameters[D];
}

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

 const size_t N = parameters.size();
 CppAD::vector<double> parametersC(N);
 for(size_t D=0; D<N;D++) {
  parametersC[D] = parameters[D];
}

Rcpp::LogicalMatrix result = hessianDenseLogical(
  parametersC, dataC, configC
  );
return(result);
}



//' @export
// [[Rcpp::export]]
Rcpp::RObject gradLogical(
  const Rcpp::NumericVector parameters,
  const Rcpp::List& data,
  const Rcpp::List& config
  ) {
  const Data   dataC(data);
  const Config configC(config);

  const size_t Nparams  = static_cast<size_t>(parameters.size());
  const size_t Nstrata  = static_cast<size_t>(dataC.Nstrata);

  // Copy parameters into a CppAD vector
  CppAD::vector<double> parametersC(Nparams);
  for (size_t j = 0; j < Nparams; ++j)
    parametersC[j] = parameters[j];

  // Output: logical pattern [parameter x strata]
  Rcpp::LogicalMatrix result(Nparams, Nstrata);

  // Optional: numeric gradients for debugging
  Rcpp::NumericMatrix resultDense;
  if (configC.debug) {
    resultDense = Rcpp::NumericMatrix(Nparams, Nstrata);
  }

  if (configC.verbose) {
    Rcpp::Rcout << "gradLogical: using Jacobian sparsity (for_jac_sparsity), single-threaded\n";
  }

  // Loop over strata (no OpenMP)
  for (size_t Dstrata = 0; Dstrata < Nstrata; ++Dstrata) {

    CppAD::vector< CppAD::AD<double> > ad_params(Nparams);
    for (size_t j = 0; j < Nparams; ++j)
      ad_params[j] = parametersC[j];

    CppAD::Independent(ad_params);

    auto latent = unpack_params(ad_params, dataC, configC);

    auto etaHere = compute_eta_for_stratum(
      Dstrata, dataC, latent.gamma, latent.beta
      );

    auto out = accumulate_contrib_for_stratum(
      Dstrata, dataC, etaHere, latent, configC
    );  // expect length-1 AD vector

    CppAD::ADFun<double> fun(ad_params, out);

    const size_t n = fun.Domain();
    const size_t m = fun.Range();    // should be 1
    if (n != Nparams || m != 1) {
      Rcpp::stop("gradLogical: fun.Domain() != Nparams or fun.Range() != 1");
    }

    CppAD::sparse_rc< CppAD::vector<size_t> > pattern_in;
    pattern_in.resize(n, n, n);   // n nonzeros
    for (size_t j = 0; j < n; ++j) {
      pattern_in.set(j, j, j);   // (row=j, col=j, index=j)
    }

    bool transpose    = false;   // want pattern_out as m × n
    bool dependency   = false;   // standard Jacobian sparsity (not dependency)
    bool internal_bool = false;  // let CppAD choose representation; updated on return

    CppAD::sparse_rc< CppAD::vector<size_t> > pattern_out;
    fun.for_jac_sparsity(
      pattern_in,
      transpose,
      dependency,
      internal_bool,
      pattern_out
      );

//    const CppAD::vector<size_t>& row = pattern_out.row();
    const CppAD::vector<size_t>& col = pattern_out.col();
    const size_t K = pattern_out.nnz();

    if (configC.verbose) {
      Rcpp::Rcout << "    sparsity nnz = " << K << "\n";
    }

    // For scalar output, row[k] should always be 0, and col[k] is the parameter index
    for (size_t k = 0; k < K; ++k) {
      size_t j = col[k];
      if (j < Nparams) {
        result(j, Dstrata) = true;
      }
    }

    // --- Optional numeric gradient for debugging ---

    if (configC.debug) {
      // Jacobian has size m * n = 1 * Nparams
      CppAD::vector<double> jac = fun.Jacobian(parametersC);
      for (size_t j = 0; j < Nparams; ++j) {
        resultDense(j, Dstrata) = jac[j];  // row = parameter, col = strata
      }
    }
  } // end strata loop

  if (configC.debug) {
    return resultDense;  // numeric gradients
  } else {
    return result;       // structural logical pattern
  }
}
