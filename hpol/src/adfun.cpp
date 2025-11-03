#include "adfun.hpp"
#include"loglikHelpers.hpp"
#include"foromp.hpp"

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

      auto latent=unpack_params(ad_params, data, config);

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


// helper: convert 1-based R indices to 0-based and keep (i >= j) lower triangle
CPPAD_TESTVECTOR( std::set<size_t> ) build_pattern_from_R(
  const Rcpp::IntegerVector& row0,
  const Rcpp::IntegerVector& col0,
  size_t n) {
 auto K = row0.size();
 CPPAD_TESTVECTOR(std::set<size_t>) pattern(n);
 
 for (size_t k = 0; k < K; ++k) {
    int ri = row0[k];   // convert to 0-based if R passed 1-based
    int cj = col0[k];
    pattern[(size_t)ri].insert((size_t)cj);
    pattern[(size_t)cj].insert((size_t)ri);
  }

  return pattern;
}


GroupPack getAdFunQ(
  const std::vector<double>& parameters,
               const Data&                data,
               const Config&              config)
{

  const Rcpp::List sparsityQ = config.sparsity["Q"];

  const size_t Nparams = parameters.size();
  const bool onlyRandom = (Nparams == data.Ngamma);

  if(config.verbose) Rcpp::Rcout << "q random " << onlyRandom << " Nparams " << Nparams;

  const Rcpp::List nsQ =  onlyRandom?sparsityQ["randomNS"]:sparsityQ["nonSymmetric"]; 
  const Rcpp::List outList = onlyRandom?sparsityQ["random"]:sparsityQ["full"];


  const Rcpp::IntegerVector rowQns = nsQ["i"]; 
  const Rcpp::IntegerVector colQns = nsQ["j"]; 
  const Rcpp::IntegerVector pNs = nsQ["p"]; 
  const Rcpp::IntegerVector matchQns = nsQ["match"]; 

  const auto pattern = build_pattern_from_R(rowQns, colQns, Nparams);

      const Rcpp::IntegerVector rowOut = outList["i"]; 
      const Rcpp::IntegerVector colOut = outList["j"]; 
      const Rcpp::IntegerVector pOut = outList["p"]; 
      const Rcpp::IntegerVector matchOut = outList["match"]; 

      std::array< std::vector<size_t>, 3 > outRowCol, nsRowCol;

      outRowCol[0] = Rcpp::as<std::vector<size_t>>(rowOut);
      outRowCol[1] = Rcpp::as<std::vector<size_t>>(colOut);
      outRowCol[2] = Rcpp::as<std::vector<size_t>>(matchOut);

      std::vector<size_t> pvec = Rcpp::as<std::vector<size_t>>(pOut);


      nsRowCol[0] = Rcpp::as<std::vector<size_t>>(rowQns);
      nsRowCol[1] = Rcpp::as<std::vector<size_t>>(colQns);
      nsRowCol[2] = Rcpp::as<std::vector<size_t>>(matchQns);

      std::vector<size_t> pvecNs = Rcpp::as<std::vector<size_t>>(pNs);


      CppAD::sparse_hessian_work work;

      CppAD::ADFun<double>  fun = adFunQ(parameters, data, config);

      GroupPack result;
      result.fun=fun;
      result.work= work;
      result.pattern=pattern;
      result.outRowCol =  outRowCol;
      result.nsRowCol = nsRowCol;
      result.outP = pvec;
      result.nsP = pvecNs;

      return(result);
}


inline std::vector<GroupPack>
getAdFun(const std::vector<double>& parameters,
                           const Data&                data,
                           const Config&              config)
{
  const Rcpp::List sparsityList = config.group_sparsity;
  const Rcpp::List strata = config.groups;
  const size_t Nparams = parameters.size();
  const bool onlyRandom = Nparams == data.Ngamma;
  const Rcpp::IntegerVector strataI = strata["i"], strataP = strata["p"];

  if(config.verbose) Rcpp::Rcout << "lik: random " << onlyRandom << " Nparams " << Nparams;

  const size_t Ngroup  = static_cast<size_t>(strataP.size() - 1);

  std::vector<GroupPack> packs(Ngroup);

  // ---- Phase 1: extract everything from R (single-thread, safe) ----
  if(sparsityList.size()) {
  for (size_t g = 0; g < Ngroup; ++g) {
    Rcpp::List sparsityHere = sparsityList[g];
    Rcpp::List second       = sparsityHere["second"];

    Rcpp::List nonSymmetric = onlyRandom ? second["randomNS"] : second["nonSymmetric"];
    Rcpp::IntegerVector Srow = nonSymmetric["i"];
    Rcpp::IntegerVector Scol = nonSymmetric["j"];
    Rcpp::IntegerVector pNs = nonSymmetric["p"];
    Rcpp::IntegerVector matchNs = nonSymmetric["match"];

    // full symmetric pattern for this group (size = Nparams)
    packs[g].pattern = build_pattern_from_R(Srow, Scol, Nparams);

    // output subset (rows/cols), copied to STL (thread-safe)
    Rcpp::List outList          = onlyRandom ? second["random"] : second["full"];
    Rcpp::IntegerVector rowOutR = outList["i"];
    Rcpp::IntegerVector colOutR = outList["j"];
    Rcpp::IntegerVector matchOutR = outList["match"];
    Rcpp::IntegerVector pOutR = outList["p"];

    packs[g].outRowCol[0] = Rcpp::as<std::vector<size_t>>(rowOutR);
    packs[g].outRowCol[1] = Rcpp::as<std::vector<size_t>>(colOutR);
    packs[g].outRowCol[2] = Rcpp::as<std::vector<size_t>>(matchOutR);
    std::vector<size_t> pvec = Rcpp::as<std::vector<size_t>>(pOutR);
    packs[g].outP = pvec;
 

    packs[g].nsRowCol[0] = Rcpp::as<std::vector<size_t>>(Srow);
    packs[g].nsRowCol[1]  = Rcpp::as<std::vector<size_t>>(Scol);
    packs[g].nsRowCol[2] = Rcpp::as<std::vector<size_t>>(matchNs);

    std::vector<size_t> pvecNs = Rcpp::as<std::vector<size_t>>(pNs);
 
    packs[g].nsP = pvecNs;

    // default-constructed work is fine; keep it per group
    packs[g].work = CppAD::sparse_hessian_work();
  }
  }

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
    const std::vector<double>& parametersC,
    const Data& dataC,
    const Config& configC)
{
  AdpackHandle h;

  if (adFun != R_NilValue) {
    // Rehydrate external pointer
    Rcpp::XPtr<std::vector<GroupPack>> xp(adFun);

    // Optional: validate class tag
    if (xp.hasAttribute("class")) {
      Rcpp::CharacterVector cls = xp.attr("class");
      if (cls.size() == 0 || std::string(cls[0]) != "adpack_ptr")
        Rcpp::stop("getAdpackFromR: external pointer class mismatch");
    }

    h.ptr = xp.get();              // we do NOT own this memory
    h.created_here = false;
  } else {
    // Build on the fly (we own it)
    auto packs = getAdFun(parametersC, dataC, configC);
    h.ptr = new std::vector<GroupPack>(std::move(packs));
    h.created_here = true;
  }

  return h;
}

