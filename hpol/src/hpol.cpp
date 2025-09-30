#include"hpol.hpp"

//#define DEBUG


#include<omp.h>


// ---- CppAD parallel hooks ----
bool in_parallel() { return omp_in_parallel() != 0; }
size_t thread_number() { return static_cast<size_t>(omp_get_thread_num()); }


using Rcpp::List; using Rcpp::S4;
using Rcpp::IntegerVector; using Rcpp::NumericVector; using Rcpp::LogicalVector;

// ----- safe getters -----
inline bool get_bool(const List& cfg, const char* key, bool def=false) {
  return cfg.containsElementNamed(key) ? Rcpp::as<bool>(cfg[key]) : def;
}
// helper to fetch int from R list with default
inline int get_int(const Rcpp::List& cfg, const char* key, int def = 1) {
  return cfg.containsElementNamed(key) ? Rcpp::as<int>(cfg[key]) : def;
}

// ----- a lightweight view of a dgCMatrix (no copies; just SEXP handles) -----
// Lightweight view for Matrix::*gCMatrix (column-compressed)
// Works for dgCMatrix (numeric x) and ngCMatrix (no x: pattern-only).
struct DgCView {
  Rcpp::IntegerVector i;    // row indices (0-based)
  Rcpp::IntegerVector p;    // column pointers (length ncol+1)
  Rcpp::NumericVector x;    // may be length 0 for ngCMatrix
  Rcpp::IntegerVector Dim;  // c(nrow, ncol)
  bool has_x;               // true if numeric/logical 'x' present

  explicit DgCView(const Rcpp::S4& obj)
  : i(obj.slot("i")),
  p(obj.slot("p")),
  Dim(obj.slot("Dim"))
  {
    // ngCMatrix has no 'x' slot; others do (dgCMatrix: numeric, lgCMatrix: logical)
    if (obj.inherits("ngCMatrix")) {
      x = Rcpp::NumericVector();   // empty
      has_x = false;
    } else {
      // coerce any existing x (numeric/logical) to NumericVector
      x = Rcpp::as<Rcpp::NumericVector>(obj.slot("x"));
      has_x = true;
    }
  }

  inline int nrow() const { return Dim[0]; }
  inline int ncol() const { return Dim[1]; }
  inline R_xlen_t nnz() const { return i.size(); } // structural nnz

  // Get numeric value for kth stored nonzero; returns 1 for pattern matrices.
  template <class T = double>
  inline T value(R_xlen_t k) const {
    return has_x ? static_cast<T>(x[k]) : static_cast<T>(1);
  }
};
// Lightweight view for Matrix::dsTMatrix (symmetric triplet), no uplo logic.
struct DgTView {
  Rcpp::IntegerVector i;    // row indices (0-based), as-is
  Rcpp::IntegerVector j;    // col indices (0-based), as-is
  Rcpp::NumericVector x;    // nonzeros, as-is
  Rcpp::IntegerVector Dim;  // c(nrow, ncol)

  explicit DgTView(const Rcpp::S4& obj)
  : i(obj.slot("i")),
  j(obj.slot("j")),
  x(obj.slot("x")),
  Dim(obj.slot("Dim"))
  {
    // (optional) assert type:
    // if (!obj.inherits("dsTMatrix")) Rcpp::stop("Expected dsTMatrix");
  }

  inline int nrow() const { return Dim[0]; }
  inline int ncol() const { return Dim[1]; }
  inline R_xlen_t nnz() const { return x.size(); }
};



// ----- all data bundled as a structure
struct Data {
  // transposed matrices from your comment
  DgTView QsansDiag;   // dgCMatrix (off-diagonals in x; diag kept separately)
  DgCView A;           // ATp
  DgCView X;           // XTp
  DgCView CC;          // cc_matrixTp

  // scalars/vectors
  NumericVector Qdiag;     // diagonals of Q
  IntegerVector map;       // theta->gamma map
  IntegerVector y;         // response

  // convenience sizes (read-only)
  const size_t Nq;         // = QsansDiag.nnz()
  const size_t Nbeta;      // = X.nrow()   (remember X is transposed)
  const size_t Ngamma;     // = A.nrow()
  const size_t Neta;       // = A.ncol()
  const size_t Nstrata;    // = CC.ncol()

  explicit Data(const List& data)
  : QsansDiag( S4(data["QsansDiag"]) )
  , A(         S4(data["ATp"]) )
  , X(         S4(data["XTp"]) )
  , CC(        S4(data["cc_matrixTp"]) )
  , Qdiag(     data["Qdiag"] )
  , map(       data["map"] )
  , y(         data["y"] )
  , Nq(        QsansDiag.nnz() )
  , Nbeta(     static_cast<size_t>(X.nrow()) )
  , Ngamma(    static_cast<size_t>(A.nrow()) )
  , Neta(      static_cast<size_t>(A.ncol()) )
  , Nstrata(   static_cast<size_t>(CC.ncol()) )
  {}
};

// ----- config bundle -----
struct Config {
  Rcpp::List list;             
  bool verbose;
  bool dirichlet;
  bool transform_theta;

  // optional numeric vectors (may be length 0)
  Rcpp::NumericVector beta;
  Rcpp::NumericVector theta;
  int num_threads;

  explicit Config(const Rcpp::List& cfg)
  : list(cfg),
  verbose(get_bool(cfg, "verbose", false)),
  dirichlet(get_bool(cfg, "dirichlet", false)),
  transform_theta(get_bool(cfg, "transform_theta", false)),
  beta(cfg.containsElementNamed("beta") ? cfg["beta"] : Rcpp::NumericVector()),
  theta(cfg.containsElementNamed("theta") ? cfg["theta"] : Rcpp::NumericVector()),
  num_threads(get_int(cfg, "num_threads", 1))           // <-- int, default 1
  {}
};

template<class Type>
struct PackedParams {
  CppAD::vector<Type> beta;     // size Nbeta
  CppAD::vector<Type> gamma;    // size Ngamma
  CppAD::vector<Type> theta;    // size Ntheta (on natural scale)
  CppAD::vector<Type> logTheta; // size Ntheta (log scale)
  size_t startGamma = 0;        // index into ad_params where gamma starts (if sourced there)
};

inline Rcpp::NumericVector get_numvec_or_warn(const Rcpp::List& cfg,
  const char* key,
  const char* what) {
  if (!cfg.containsElementNamed(key)) {
    Rcpp::warning("%s missing from config", what);
    return Rcpp::NumericVector();
  }
  return Rcpp::as<Rcpp::NumericVector>(cfg[key]);
}

template<class Type>
inline void fill_from_numvec(CppAD::vector<Type>& dst,
 const Rcpp::NumericVector& src,
 const char* what) {
  if (src.size() == 0) return; // was missing; keep dst as initialized
  if (static_cast<size_t>(src.size()) != dst.size()) {
    Rcpp::stop("%s has wrong length: expected %zu, got %d",
     what, dst.size(), src.size());
  }
  for (size_t i = 0; i < dst.size(); ++i) dst[i] = Type(src[i]);
}


template<class Type>
PackedParams<Type>
unpack_params(const CppAD::vector<Type>& ad_params,
  const Data& data,
  const Config& cfg)
{
  using CppAD::exp;
  using CppAD::log;

  const size_t Nbeta  = data.Nbeta;
  const size_t Ngamma = data.Ngamma;
  const size_t Ntheta_base =
  (data.map.size() > 0 ? static_cast<size_t>(Rcpp::max(data.map)) + 1 : 0u);
  const size_t Ntheta = Ntheta_base + (cfg.dirichlet ? 1u : 0u);

  PackedParams<Type> out;
  out.beta     = CppAD::vector<Type>(Nbeta);
  out.gamma    = CppAD::vector<Type>(Ngamma);
  out.theta    = CppAD::vector<Type>(Ntheta);
  out.logTheta = CppAD::vector<Type>(Ntheta);
  out.startGamma = 0;

  const size_t n = ad_params.size();

  if (n == Ngamma) {
    // gamma only in ad_params; beta and theta must come from cfg
    for (size_t i = 0; i < Ngamma; ++i)
      out.gamma[i] = ad_params[i];

    // beta from cfg (must be present and correct length)
    if (static_cast<size_t>(cfg.beta.size()) != Nbeta)
      Rcpp::stop("beta has wrong length: expected %zu, got %d",
       Nbeta, cfg.beta.size());
    for (size_t i = 0; i < Nbeta; ++i)
      out.beta[i] = static_cast<Type>(cfg.beta[i]);

    // theta/logTheta from cfg (must be present and correct length)
    if (static_cast<size_t>(cfg.theta.size()) != Ntheta)
      Rcpp::stop("theta has wrong length: expected %zu, got %d",
       Ntheta, cfg.theta.size());

    if (cfg.transform_theta) {
      // cfg.theta contains log-theta
      for (size_t i = 0; i < Ntheta; ++i) {
        out.logTheta[i] = static_cast<Type>(cfg.theta[i]);
        out.theta[i]    = exp(out.logTheta[i]);
      }
    } else {
      // cfg.theta contains natural theta
      for (size_t i = 0; i < Ntheta; ++i) {
        out.theta[i]    = static_cast<Type>(cfg.theta[i]);
        out.logTheta[i] = log(out.theta[i]);
      }
    }
  } else {
    // Everything in ad_params: [beta | gamma | theta-or-logtheta]
    const size_t expected = Nbeta + Ngamma + Ntheta;
    if (n != expected) {
      Rcpp::stop("parameters has wrong size: expected %zu (Nbeta+Ngamma+Ntheta), got %zu",
       expected, n);
    }

    // beta
    for (size_t i = 0; i < Nbeta; ++i)
      out.beta[i] = ad_params[i];

    // gamma
    out.startGamma = Nbeta;
    for (size_t i = 0; i < Ngamma; ++i)
      out.gamma[i] = ad_params[out.startGamma + i];

    // theta or log-theta
    const size_t startTheta = Nbeta + Ngamma;
    if (cfg.transform_theta) {
      // ad_params contains log-theta
      for (size_t i = 0; i < Ntheta; ++i) {
        out.logTheta[i] = ad_params[startTheta + i];
        out.theta[i]    = exp(out.logTheta[i]);
      }
    } else {
      // ad_params contains natural theta
      for (size_t i = 0; i < Ntheta; ++i) {
        out.theta[i]    = ad_params[startTheta + i];
        out.logTheta[i] = log(out.theta[i]);
      }
    }
  }

  return out;
}




template<class Type>
CppAD::vector<Type>  objectiveFunctionInternal(
  CppAD::vector<Type> params_input, 
  Rcpp::List dataList, 
  Rcpp::List configList 
  ) {


  Data   data(dataList);
  Config cfg(configList);

  auto latent = unpack_params<Type>(params_input, data, cfg);
  // elements beta, gamma, theta, logtheta


  CppAD::vector<Type> eta(data.Neta, Type(0));
  CppAD::vector<Type> etaLogSum(data.Nstrata);
  CppAD::vector<Type> gammaScaled(data.Ngamma);

  CppAD::vector<Type> minusLogDens(1,0);


  Type nu = latent.theta[latent.theta.size()-1],
  logSqrtNu = nu/ 2,
  oneOverSqrtNu = exp(-logSqrtNu),
  lgammaOneOverSqrtNu = lgamma_ad(oneOverSqrtNu);


#ifdef DEBUG
  Rcpp::Rcout << "nu " << nu << " logSqrtNU " << logSqrtNu << 
  " oneOverSqrtNu " << oneOverSqrtNu << " lgammaOneOverSqrtNu " <<
  lgammaOneOverSqrtNu << std::endl << "theta ";
  for(int Dtheta =0;Dtheta < latent.theta.size();++Dtheta) {
    Rcpp::Rcout << Dtheta << " " << latent.theta[Dtheta] << " " << latent.logTheta[Dtheta] << std::endl;
  }
#endif    

  size_t num_threads;
  if constexpr (std::is_same<Type, double>::value) {
        // ordinary numeric run → use requested number of threads
    num_threads=cfg.num_threads;
  } else if constexpr (std::is_same<Type, CppAD::AD<double>>::value) {
        // AD recording → force serial execution
    num_threads =1;
  }
  omp_set_num_threads(num_threads);


  CppAD::vector<Type> randomContributionDiag(num_threads);
  CppAD::vector<Type>  offdiagQ(num_threads);
  CppAD::vector<Type>  loglik(num_threads);


      #pragma omp parallel
  { 
    const int tid=omp_get_thread_num();

    randomContributionDiag[tid] = Type(0);
    offdiagQ[tid]               = Type(0);
    loglik[tid]                 = Type(0);

  // eta = X * beta + A * gamma   
    #pragma omp for nowait
    for (size_t Deta = 0; Deta < data.Neta; ++Deta) {
      const int p0x = data.X.p[Deta];
      const int p1x = data.X.p[Deta + 1];
      const int p0a = data.A.p[Deta];
      const int p1a = data.A.p[Deta + 1];    
      eta[Deta]=Type(0);

  // X contribution: columns are eta indices  
      for (int k = p0x; k < p1x; ++k) {
      const int    r = data.X.i[k];          // beta row index
      const Type   v = Type(data.X.x[k]);    // value
      eta[Deta]   += v * latent.beta[r];
    }

  // A contribution
    for (int k = p0a; k < p1a; ++k) {
      const int    r = data.A.i[k];          // gamma row index
      const Type   v = Type(data.A.x[k]);    // value
      eta[Deta]   += v * latent.gamma[r];
    }
  }

  // log(|Q|) + 0.5 * gamma^T Q gamma

  
  // diagonals
  #pragma omp for 
  for(size_t D=0;D<data.Ngamma;D++) {
    size_t mapHere = data.map[D];

    gammaScaled[D] = latent.gamma[D] / latent.theta[mapHere];
    randomContributionDiag[tid] += latent.logTheta[mapHere] +
    0.5*gammaScaled[D]*gammaScaled[D]*data.Qdiag[D];
  }

// need to make sure all gammaScaled and etas are created before moving on
#pragma omp barrier


#ifdef DEBUG
  Rcpp::Rcout << "Q offdiag " << data.Nq << std::endl;
#endif    

    // Q offdiag    
      #pragma omp for nowait
  for(size_t D = 0; D < data.Nq; D++) {
    offdiagQ[tid] += gammaScaled[data.QsansDiag.i[D]] * gammaScaled[data.QsansDiag.j[D]] * data.QsansDiag.x[D];
  }


#ifdef DEBUG
  Rcpp::Rcout << "sum log etas\n";
#endif
  // calculate log(sum(exp(eta_i))) within strata

    #pragma omp for nowait
  for (size_t i = 0; i < data.Nstrata; i++) {
    size_t startHere = data.CC.p[i], Nhere = data.CC.p[i+1];
    size_t NinStrata = Nhere - startHere;
    CppAD::vector<Type> etaHere(NinStrata);


    // loop through row i of CCmatrix
    // this is j=startHere
    for(size_t j0=0, j=startHere; j < Nhere; j++,j0++) {
      etaHere[j0] = eta[data.CC.i[j]];
    }
#ifdef USEATOMICS    
    // use the atomic formula stuff with analytical derivatives
    etaLogSum[i] = logspace_add_n(etaHere);
#else 
    size_t max_idx = 0;
    for(size_t Didx = 1; Didx < etaHere.size(); ++Didx) {
      if(etaHere[Didx] > etaHere[max_idx]) {
        max_idx = Didx;
      }
    }
//    double max_value = CppAD::Value(etaHere[max_idx]);
    double max_value;
    if constexpr (std::is_same<Type, double>::value) {
      max_value = etaHere[max_idx];
    }  else {
      max_value = CppAD::Value(etaHere[max_idx]);
    } 

    Type sumexp=0.0;
    for(size_t j = 0; j < etaHere.size(); ++j) {
      sumexp += exp(etaHere[j] - max_value);
    }
    etaLogSum[i] = max_value + log(sumexp);
#endif    




// for data contribution

    Type  contrib = 0.0;
    int sumY = 0;

    for(size_t j=startHere; j < Nhere; j++) {
      size_t idx = data.CC.i[j];
      sumY += data.y[idx];     
      Type etaMinusLogSumMu = eta[idx] - etaLogSum[i];
      Type  muBarDivSqrtNu = exp(etaMinusLogSumMu - logSqrtNu);
#ifdef TEST
      contrib -= (data.y[idx] - eta[idx])*(data.y[idx] - eta[idx])/(nu*nu);
#else      
      if(cfg.dirichlet) {
        contrib += lgamma_ad(data.y[idx] + muBarDivSqrtNu) - 
        lgamma_ad(muBarDivSqrtNu);
      } else {
        contrib += data.y[idx] * etaMinusLogSumMu;
      }
#endif
      
#ifdef EVALCONSTANTS
      contrib -= lgamma(data.y[idx] + 1);
#endif
    } // j obs in strata

    if(cfg.dirichlet) {
      contrib += lgammaOneOverSqrtNu - lgamma_ad(oneOverSqrtNu + sumY);
    }
#ifdef EVALCONSTANTS
    contrib += lgamma(sumY + 1);
#endif

    loglik[tid] += contrib;
  } // loop through strata

} //parellel block

Type randomContributionDiagS = Type(0);
Type  offdiagQS= Type(0);
Type  loglikS = Type(0);



for (int t = 0; t < num_threads; ++t) {
#ifdef DEBUG
  Rcpp::Rcout << "thread " << t << "logLik " << loglik[t] << " offdiag " << offdiagQ[t] <<
  " diag " << randomContributionDiag[t] << std::endl;
#endif 
  randomContributionDiagS += randomContributionDiag[t];
  offdiagQS += offdiagQ[t];
  loglikS += loglik[t];
}

#ifdef DEBUG
Rcpp::Rcout << "logLik " << loglikS << " offdiag " << offdiagQS <<
" diag " << randomContributionDiagS << std::endl;
#endif 


minusLogDens[0] =  - loglikS + offdiagQS  + randomContributionDiagS;
//  etaLogSum[0]  + etaLogSum[1]


//    Rcpp::Rcout << " " << etaLogSum[0] << " " << etaLogSum[1] << "\n";

#ifdef EVALCONSTANTS
minusLogDens[0] += Ngamma * HALFLOGTWOPI;
//  minusLogDens -= 0.5*logdet(Q);
if (config.containsElementNamed("halfLogDetQ")) {
  minusLogDens[0] -= Rcpp::as<double>(config["halfLogDetQ"]);
}
#endif

#ifdef DEBUG
Rcpp::Rcout << "all done " << minusLogDens[0] << std::endl;
#endif 

return minusLogDens;
}


/* Parameters: 
 X A cc_matrix: dgRMatrix, cc_matrix can be lgRMatrix or ngRMatrix (with MatrixExtra package)
 QsansDiag, Qdiag dsTMatrix and vector for precision matrix
 map, y integer vectors 
 config: hesMax integer, maximum number of non-zero hessian elements
 */

//' @export
// [[Rcpp::export]]
double objectiveFunctionNoDiff(
  Rcpp::NumericVector parameters, 
  Rcpp::List data, 
  Rcpp::List config
  ) {

  size_t Nparams = parameters.size();
  CppAD::vector<double> parametersC(Nparams);
  for (size_t D = 0; D < Nparams; D++) {
    parametersC[D] = parameters[D];  // Initialize CppAD variables
  }
  auto result = objectiveFunctionInternal<double>(parametersC, data, config);
  return result[0];
}

//' @export
// [[Rcpp::export]]
Rcpp::List objectiveFunctionC(
  Rcpp::NumericVector parameters, 
  Rcpp::List data, 
  Rcpp::List config
  ) {

  bool verbose = false;
  if (config.containsElementNamed("verbose")) {
    verbose = Rcpp::as<bool>(config["verbose"]);
  }
  bool dense = false;
  if (config.containsElementNamed("dense")) {
    dense = Rcpp::as<bool>(config["dense"]);
  }
  size_t maxDeriv = 2;
  if (config.containsElementNamed("maxDeriv")) {
    maxDeriv = Rcpp::as<int>(config["maxDeriv"]);
  }
  int num_threads = 1;
  if (config.containsElementNamed("num_threads"))
    num_threads = Rcpp::as<int>(config["num_threads"]);


  size_t Nparams = parameters.size();
  int NparamsI = static_cast<int>(Nparams);
  int hesMax = NparamsI*NparamsI/2;
  if (config.containsElementNamed("hesMax")) {
    hesMax = Rcpp::as<int>(config["hesMax"]);
  }



  CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
  for (size_t D = 0; D < Nparams; D++) {
    ad_params[D] = parameters[D];  // Initialize CppAD variables
  }
  CppAD::Independent(ad_params);  // Tell CppAD these are inputs for differentiation
  
  if (verbose ) {
    Rcpp::Rcout << "eval";
  }
  CppAD::vector<CppAD::AD<double>> y = objectiveFunctionInternal(ad_params, data, config);
  if (verbose ) {
    Rcpp::Rcout << "y " << y[0] << "\n";
  }
  CppAD::ADFun<double> fun(ad_params, y);

  std::vector<double> x_val(Nparams);
  for (size_t i = 0; i < Nparams; ++i) { 
    x_val[i] = parameters[i];
  }

  if (verbose ) {
    Rcpp::Rcout << "forward: ";
  }

  std::vector<double> y_val(1);
  y_val = fun.Forward(0, x_val);

  if (verbose ) {
    Rcpp::Rcout << "f0 " << y_val[0] << "\n";
  }


  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("value") = y_val[0]
    );

  if(maxDeriv == 0) {
    return result;
  }

  if (verbose ) {
    Rcpp::Rcout << "grad ";
  }
  std::vector<double> grad = fun.Jacobian(x_val);
  if (verbose ) {
    Rcpp::Rcout << ".\n";
  }
  result["grad"] = grad;
  if(maxDeriv == 1) {
    return result;
  }

// hessian
  if (verbose ) {
    Rcpp::Rcout << "hess " << num_threads << " threads\n";
  }

  CppAD::vector< std::set<size_t> >  sparsity(Nparams);
  Rcpp::IntegerVector Hrow, Hcol, Hp;
  Rcpp::NumericVector Hvalue;
  bool haveSparsity;

  Rcpp::IntegerVector Hstart(num_threads,0), Hend(num_threads,0);

  Rcpp::NumericMatrix denseHessianOut;
  if(dense) {
    denseHessianOut = Rcpp::NumericMatrix(Nparams, Nparams);
  }

    // if have sparsity
  if (config.containsElementNamed("sparsity")) {
    if (verbose ) {
      Rcpp::Rcout << "given sparsity pattern\n";
    }      
    haveSparsity=TRUE;

    Rcpp::List sparsityR, sparsityR2, sparsityR1 = config["sparsity"];
    if (sparsityR1.containsElementNamed("second")) {
        // list with second and third sparsity
      sparsityR2 = sparsityR1["second"];
      if (config.containsElementNamed("theta")) {
        if (verbose ) {
          Rcpp::Rcout << "sparsity only for random effects\n";
        }   
        sparsityR = sparsityR2["random"];
      } else {
        if (verbose ) {
          Rcpp::Rcout << "full parameter sparsity\n";
        }   
        sparsityR = sparsityR2["full"];
      }
    } else {
      if (verbose ) {
        Rcpp::Rcout << "sparsity i j p provided\n";
      }   
      sparsityR = sparsityR1;
    }

    Hrow = sparsityR["i"]; 
    Hcol = sparsityR["j"];
    Hp = sparsityR["p"];
    Hvalue = Rcpp::NumericVector(Hrow.size());


    // where each thread starts and ends in hessian
    int NperThread = (Hvalue.size() + num_threads - 1) / num_threads;

    if (verbose ) {
      Rcpp::Rcout << "NperThread " << NperThread << 
      " non-zeros " << Hvalue.size() << "\n";
    }   


  // start and end points for each thread
    Hstart = NperThread * Rcpp::seq(0, num_threads - 1);
  Hstart = Rcpp::pmin(Hstart, Hrow.size());  // clamp
  Hend = Rcpp::pmin(Hstart + NperThread, Hrow.size());
  



      // convert to vector set
// Loop over columns (variables)
  for(int D=0;D<Hrow.size();++D) {
    int Drow = Hrow[D], Dcol=Hcol[D];
    sparsity[Drow].insert(Dcol);    
    if(Drow != Dcol) sparsity[Dcol].insert(Drow);   
  }


    } else { // if don't have sparsity
      if (verbose ) {
        Rcpp::Rcout << "compute dense hessian hesMax " << hesMax << "\n";
      }
      haveSparsity = FALSE;
      Hrow = Rcpp::IntegerVector(hesMax, -1); 
      Hcol = Rcpp::IntegerVector(hesMax, -1); 
      Hvalue = Rcpp::NumericVector(hesMax);
      int NperThread = Hvalue.size()/num_threads;
      Hend = Hstart = NperThread * Rcpp::seq(0, num_threads - 1); 
    }



    omp_set_num_threads(num_threads);
    CppAD::thread_alloc::parallel_setup(num_threads, in_parallel, thread_number);
    CppAD::thread_alloc::hold_memory(true);

    // Replicate fun object for each thread
    std::vector<CppAD::ADFun<double>> fun_threads(num_threads);
    for (int i = 0; i < num_threads; ++i) {
      fun_threads[i] = fun;
    }

    const double eps = 1e-12;
//    int hindex = 0, hindex_thread=0;



    #pragma omp parallel
    {


      const int tid = thread_number();
      std::vector<double> w(1, 1.0);

      // 
      if(haveSparsity) {
        // rounding up
        int Nelements = Hend[tid] - Hstart[tid];
        std::vector<double> thread_Hvalue(Nelements);
        std::vector<int> thread_Hrow(Nelements), thread_Hcol(Nelements);

        for(int D=0,Dorig=Hstart[tid];D<Nelements;D++,Dorig++) {
          thread_Hrow[D] = Hrow[Dorig];
          thread_Hcol[D] = Hcol[Dorig];
        }

#ifdef DEBUG
        Rcpp::Rcout << "t" << tid << "N" << Nelements <<
        "s" << Hstart[tid] << "e" << Hend[tid] <<
        "\n";
#endif 

        CppAD::sparse_hessian_work work; 

        fun_threads[tid].SparseHessian(
          x_val, w, sparsity, 
          thread_Hrow, thread_Hcol, thread_Hvalue, 
          work);

// copy elements to result vector
        for(
          int DfromZero=0,DfromStart=Hstart[tid];
          DfromZero < thread_Hvalue.size(); //DfromStart < Hend[tid]; 
          ++DfromZero,++DfromStart
          ) {
          Hvalue[DfromStart] = thread_Hvalue[DfromZero];
        if(dense) {
          denseHessianOut(thread_Hrow[DfromZero], 
            thread_Hcol[DfromZero]) = 
          denseHessianOut(thread_Hcol[DfromZero], 
            thread_Hrow[DfromZero]) = thread_Hvalue[DfromZero];
        }

      }

      } else { // no sparsity
        // eval dense hessian, save non-zeros


#ifdef DEBUG
        Rcpp::Rcout << "t" << tid << "\n";
#endif 

        std::vector<double> u(Nparams, 0.0);

    #pragma omp for
        for (int j = 0; j < NparamsI; j++) {
          std::fill(u.begin(), u.end(), 0.0);
          u[j] = 1.0;
          fun_threads[tid].Forward(0, x_val);
          fun_threads[tid].Forward(1, u);
          std::vector<double> ddw = fun_threads[tid].Reverse(2, w);

          for (int irow = j; irow < NparamsI; ++irow) {
            double dhere = ddw[2 * irow + 1];
            if(dense) {
              denseHessianOut(irow, j) = dhere;
            }
            if (!CppAD::NearEqual(dhere, 0.0, eps, eps)) {
              Hvalue[Hend[tid]] = dhere;
              Hrow[Hend[tid]] = irow;
              Hcol[Hend[tid]] = j;
              if(Hend[tid]<hesMax) Hend[tid]++;
            } // if dhere not zero
          } // irow
      } // j column
    } // no sparsity pattern

  } // parallel


#ifdef UNDEF
  CppAD::thread_alloc::parallel_setup(1, nullptr, nullptr);
  CppAD::thread_alloc::hold_memory(false);
  for (int i = 1; i < num_threads; ++i)
    CppAD::thread_alloc::free_available(i);
  CppAD::thread_alloc::free_available(0);
#endif

/*Rcpp::List hessianR = Rcpp::List::create(
  Rcpp::Named("i") = Hrow,
  Rcpp::Named("j") = Hcol,
  Rcpp::Named("x") = Hvalue,
  Rcpp::Named("start") = Hstart,
  Rcpp::Named("end") = Hend
);*/

  Rcpp::S4 hessianR = make_TMatrix(
    Hvalue, Hrow, Hcol, 
    Nparams); 

  result["hessian"] = hessianR;

if(dense) {
  result["denseHessian"] = denseHessianOut;
}


  if (verbose ) { 
    Rcpp::Rcout << " done." << std::endl;
  }

  return result;
}

