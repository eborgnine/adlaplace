#ifndef CONFIGOBJ_HPP
#define CONFIGOBJ_HPP

#include"adlaplace/matrixUtils.hpp"

// ----- config bundle -----
struct Config {
  bool verbose;
  bool transform_theta;
  int num_threads=1;
  std::vector<std::vector<int>> hessianIP = std::vector<std::vector<int>>(2);
  std::vector<std::vector<int>> hessianIPLower = std::vector<std::vector<int>>(3);


  std::vector<double> beta;

  // always available:
  std::vector<double> theta;     // natural scale  (theta > 0)
  std::vector<double> logTheta;  // log(theta)
  double nbSize;
  double logNbSize;
  double lgammaNbSize;
  double sizeLogSize;


  explicit Config(const Rcpp::List& cfg)
  : verbose(get_bool(cfg, "verbose", false)),
  transform_theta(get_bool(cfg, "transform_theta", false)),
  num_threads(get_int(cfg, "num_threads", 1)),           
  beta(get_numvec_copy(cfg, "beta")) 
  {

    if(cfg.containsElementNamed("hessian")) {
  const Rcpp::S4 sm = cfg["hessian"];
  const Rcpp::IntegerVector i = sm.slot("i");
  const Rcpp::IntegerVector p = sm.slot("p");
  hessianIP[0] = std::vector<int>(i.begin(), i.end());
  hessianIP[1] = std::vector<int>(p.begin(), p.end());
    }
    if(cfg.containsElementNamed("hessianL")) {
  const Rcpp::S4 sm = cfg["hessianL"];
  const Rcpp::IntegerVector i = sm.slot("i");
  const Rcpp::IntegerVector p = sm.slot("p");
  const Rcpp::IntegerVector x = sm.slot("x");
  hessianIPLower[0] = std::vector<int>(i.begin(), i.end());
  hessianIPLower[1] = std::vector<int>(p.begin(), p.end());
  hessianIPLower[2] = std::vector<int>(x.begin(), x.end());
    }


    std::vector<double> theta_input = get_numvec_copy(cfg, "theta");

    if (transform_theta) {
      // User supplied log(theta); store as logTheta, derive theta
      logTheta = theta_input;
      theta.resize(logTheta.size());
      for (std::size_t i = 0; i < logTheta.size(); ++i) {
        theta[i] = std::exp(logTheta[i]);
      }
    } else {
      // User supplied theta on natural scale; store as theta, derive logTheta
      theta = theta_input;
      logTheta.resize(theta.size());
      for (std::size_t i = 0; i < theta.size(); ++i) {
        // you may want to guard theta[i] > 0 here
        logTheta[i] = std::log(theta[i]);
      }
    }
    // last element is 1/sqrt(nbSize)
    logNbSize = -0.5*logTheta[logTheta.size()-1];
    nbSize = std::exp(logNbSize);
    lgammaNbSize = std::lgamma(nbSize);    
    sizeLogSize = nbSize * logNbSize;
    if(verbose) {
      Rcpp::Rcout << "last theta " << theta[theta.size()-1] <<
      " log theta " << logTheta[theta.size()-1] <<
      " nbSize " <<
      nbSize << " lgammaNbSize " << lgammaNbSize << "\n";
    }
  }
};
#endif