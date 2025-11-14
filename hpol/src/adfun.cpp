#include "adfun.hpp"
#include"loglikHelpers.hpp"
#include"foromp.hpp"



/* likelihood part */
CppAD::ADFun<double> adFunGroup(
      const CppAD::vector<double> & parameters,  
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

      auto latent=unpack_params(ad_params, data, config);

//      Rcpp::Rcout << "ad beta " << latent.beta[0] << " gamma " << latent.gamma[0] << "\n";

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
  const CppAD::vector<double> & parameters,  
  const Data& data,
  const Config& config) {

  const size_t Nparams = parameters.size();
  const size_t Nq = data.Qdiag.size();

  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
  for(size_t D=0;D<Nparams;D++) {
    ad_params[D] = parameters[D];
  }
  CppAD::Independent(ad_params);   

  auto latent=unpack_params<CppAD::AD<double>>(ad_params, data, config);

  CppAD::vector<CppAD::AD<double>> gammaScaled(data.Ngamma);
  CppAD::vector<CppAD::AD<double>> result(1);

  CppAD::AD<double> result0=0;
  const Rcpp::IntegerVector map = data.map;
  for (size_t D = 0; D < Nq; ++D) {
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


inline CppAD::sparse_rc< CppAD::vector<size_t> > build_hessian_pattern_from_pairs(
  const Rcpp::IntegerVector& row,
  const Rcpp::IntegerVector& col,
  size_t n) {
  const size_t K = row.size();
  CppAD::sparse_rc< CppAD::vector<size_t> > pat;
  pat.resize(n, n, K);
  bool warnRC=true;
  for (size_t k = 0; k < K; ++k) {
    const int r = row[k], c=col[k];
    if(r > c && warnRC) {
      warnRC = false;
      Rcpp::Rcout << "entry lower, should be uppper triangle\n";
    }
    pat.set(k, r, c);
  }
  return pat;
}


inline CppAD::sparse_rc< CppAD::vector<size_t> > build_gradient_pattern_from_R(
  const Rcpp::IntegerVector& grad_index, const size_t N)
{
  const size_t K = static_cast<size_t>(grad_index.size());
  std::set<size_t> idx;

  for (size_t k = 0; k < K; ++k) {
    size_t c = static_cast<size_t>(grad_index[k]);
    idx.insert(c);
  }


  CppAD::sparse_rc< CppAD::vector<size_t> > pattern;
  pattern.resize(1, N, idx.size());
  size_t t = 0;
  for (size_t j : idx) pattern.set(t++, 0, j);

  return pattern;
}


GroupPack getAdFunQ(
  const CppAD::vector<double>& parameters,
  const Data&                data,
  const Config&              config)
{

  const Rcpp::List sparsityQ = config.sparsity["Q"];

  const size_t Nparams = parameters.size();
  const bool onlyRandom = (Nparams == data.Ngamma);

//  if(config.verbose) Rcpp::Rcout << "q random " << onlyRandom << " Nparams " << Nparams;

  const Rcpp::List outList = onlyRandom?sparsityQ["random"]:sparsityQ["full"];
  const Rcpp::IntegerVector rowOut = outList["i"]; 
  const Rcpp::IntegerVector colOut = outList["j"]; 

  CppAD::ADFun<double>  fun = adFunQ(parameters, data, config);

  GroupPack result;
  result.fun=fun;
  result.work_hess= CppAD::sparse_hes_work();
  result.pattern_hess=build_hessian_pattern_from_pairs(rowOut, colOut, Nparams);
  result.out_hess = CppAD::sparse_rcv< CppAD::vector<size_t>, CppAD::vector<double> >(result.pattern_hess);

  const Rcpp::List nsQ =  onlyRandom?sparsityQ["randomNS"]:sparsityQ["nonSymmetric"]; 
  const Rcpp::IntegerVector rowQns = nsQ["i"]; 
  const Rcpp::IntegerVector colQns = nsQ["j"]; 
  const Rcpp::IntegerVector pNs = nsQ["p"]; 
  const Rcpp::IntegerVector matchQns = nsQ["match"]; 
  const Rcpp::IntegerVector pOut = outList["p"]; 
  const Rcpp::IntegerVector matchOut = outList["match"]; 

  std::vector<size_t> pvec = Rcpp::as<std::vector<size_t>>(pOut);
  std::vector<size_t> pvecNs = Rcpp::as<std::vector<size_t>>(pNs);

  std::array< std::vector<size_t>, 3 > outRowCol, nsRowCol;

  outRowCol[0] = Rcpp::as<std::vector<size_t>>(rowOut);
  outRowCol[1] = Rcpp::as<std::vector<size_t>>(colOut);
  outRowCol[2] = Rcpp::as<std::vector<size_t>>(matchOut);

  nsRowCol[0] = Rcpp::as<std::vector<size_t>>(rowQns);
  nsRowCol[1] = Rcpp::as<std::vector<size_t>>(colQns);
  nsRowCol[2] = Rcpp::as<std::vector<size_t>>(matchQns);

  result.outRowCol =  outRowCol;
  result.nsRowCol = nsRowCol;
  result.outP = pvec;
  result.nsP = pvecNs;


  return(result);
}



inline std::vector<GroupPack>
getAdFun(
  const CppAD::vector<double>& parameters,
  const Data&                data,
  const Config&              config)
{
  const Rcpp::List sparsityList = config.group_sparsity;
  const Rcpp::List strata = config.groups;
  const size_t Nparams = parameters.size();
  const bool onlyRandom = Nparams == data.Ngamma;
  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];


  const size_t Ngroup  = static_cast<size_t>(strataP.size() - 1);

  if(config.verbose) {
    Rcpp::Rcout << "lik: random " << onlyRandom << " Nparams " << Nparams << " Ngroup " << Ngroup <<
       "\n";
  }

  std::vector<GroupPack> packs(Ngroup);

  if(sparsityList.size()) {
  for (size_t g = 0; g < Ngroup; ++g) {
    const Rcpp::List sparsityHere = sparsityList[g];

    // first
    const bool haveFirst = sparsityHere.containsElementNamed("first");
    const Rcpp::List first = haveFirst?sparsityHere["first"]:Rcpp::List();
    const Rcpp::IntegerVector sparseGrad = haveFirst?(onlyRandom?first["random"]:first["full"]):Rcpp::IntegerVector();

    if(haveFirst) {
      packs[g].work_grad = CppAD::sparse_jac_work();
      packs[g].pattern_grad = build_gradient_pattern_from_R(sparseGrad, Nparams);
      packs[g].out_grad = CppAD::sparse_rcv< 
          CppAD::vector<size_t>, CppAD::vector<double> 
        >(packs[g].pattern_grad); 
    }

    // second
    const Rcpp::List second       = sparsityHere["second"];
    const Rcpp::List outList          = onlyRandom ? second["random"] : second["full"];
    const Rcpp::List nonSymmetric = onlyRandom ? second["randomNS"] : second["nonSymmetric"];
    

    const Rcpp::IntegerVector rowOutR = outList["i"];
    packs[g].outRowCol[0] = Rcpp::as<std::vector<size_t>>(rowOutR);
    const Rcpp::IntegerVector Srow = nonSymmetric["i"];
    packs[g].nsRowCol[0] = Rcpp::as<std::vector<size_t>>(Srow);
    
    const Rcpp::IntegerVector colOutR = outList["j"];
    packs[g].outRowCol[1] = Rcpp::as<std::vector<size_t>>(colOutR);
    const Rcpp::IntegerVector Scol = nonSymmetric["j"];
    packs[g].nsRowCol[1]  = Rcpp::as<std::vector<size_t>>(Scol);

    const Rcpp::IntegerVector matchOutR = outList["match"];
    packs[g].outRowCol[2] = Rcpp::as<std::vector<size_t>>(matchOutR);
    const Rcpp::IntegerVector matchNs = nonSymmetric["match"];
    packs[g].nsRowCol[2] = Rcpp::as<std::vector<size_t>>(matchNs);

    const Rcpp::IntegerVector pOutR = outList["p"];
    packs[g].outP = Rcpp::as<std::vector<size_t>>(pOutR);
    const Rcpp::IntegerVector pNs = nonSymmetric["p"];
    packs[g].nsP = Rcpp::as<std::vector<size_t>>(pNs);


    packs[g].work_hess = CppAD::sparse_hes_work();
    packs[g].pattern_hess = build_hessian_pattern_from_pairs(rowOutR, colOutR, Nparams);
    packs[g].out_hess = CppAD::sparse_rcv< 
        CppAD::vector<size_t>, CppAD::vector<double> 
      >(packs[g].pattern_hess);

  } //gropu loop
  } // sparsity list size



  // ---- Phase 2: build ADFun per group (parallel, no Rcpp touched) ----
  omp_set_num_threads(config.num_threads);
  CppAD::thread_alloc::parallel_setup(
    config.num_threads,
    [](){ return in_parallel_wrapper(); },
    [](){ return static_cast<size_t>(thread_num_wrapper()); }
  );

#pragma omp parallel for
  for (size_t g = 0; g < Ngroup; ++g) {
    // Create the taped function for this group's strata range
    CppAD::ADFun<double> f =
      adFunGroup(parameters, data, config, strataI,
                 static_cast<size_t>(strataP[g]),
                 static_cast<size_t>(strataP[g + 1]));
    packs[g].fun = std::move(f); // move-assign into the bundle
  }

  return packs;
}




// ---- helper to rehydrate or create adpack (no Nullable<SEXP>) ----
AdpackHandle getAdpackFromR(
    SEXP adFun,                                   // pass R_NilValue if none
    const CppAD::vector<double>& parametersC,
    const Data& dataC,
    const Config& configC)
{
  AdpackHandle h;

  if (adFun != R_NilValue) {
    // Rehydrate external pointer
    Rcpp::XPtr<std::vector<GroupPack>> xp(adFun);

    // Optional: validate class tag
    if (xp.hasAttribute("class")) {
      if(configC.verbose) {
        Rcpp::Rcout << "existing functions\n";
      }
      Rcpp::CharacterVector cls = xp.attr("class");
      if (cls.size() == 0 || std::string(cls[0]) != "adpack_ptr")
        Rcpp::stop("getAdpackFromR: external pointer class mismatch");
    }

    h.ptr = xp.get();              // we do NOT own this memory
    h.created_here = false;
  } else {
      if(configC.verbose) {
        Rcpp::Rcout << "creating new functions\n";
      }
    // Build on the fly (we own it)
    auto packs = getAdFun(parametersC, dataC, configC);
    h.ptr = new std::vector<GroupPack>(std::move(packs));
    h.created_here = true;
  }

  return h;
}

