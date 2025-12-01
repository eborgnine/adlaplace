#ifndef CONFIGOBJ_HPP
#define CONFIGOBJ_HPP

#include"adlaplace/matrixUtils.hpp"
#include <cppad/cppad.hpp>




// ----- config bundle -----
struct Config {
  bool verbose;
  bool transform_theta;
  int num_threads=1;

  CppAD::vector<double> beta;
  CppAD::vector<double> start_gamma;
  CppAD::vector<double> theta;    


  std::vector<std::vector<int>> hessianIP = std::vector<std::vector<int>>(2);
  std::vector<std::vector<int>> hessianIPLower = std::vector<std::vector<int>>(3);

  Rcpp::List group_sparsity;


  explicit Config(const Rcpp::List& cfg)
  : verbose(get_bool(cfg, "verbose", false)),
  transform_theta(get_bool(cfg, "transform_theta", false)),
  num_threads(get_int(cfg, "num_threads", 1))     
  {

    const Rcpp::NumericVector start_gamma_here = cfg["start_gamma"];
    const size_t Nparams = start_gamma_here.size();
    start_gamma = CppAD::vector<double>(Nparams);
    for(size_t D=0;D<Nparams;D++) {
      start_gamma[D] = start_gamma_here[D];
    }

    std::vector<double> theta_input = get_numvec_copy(cfg, "theta");
    theta = CppAD::vector<double>(theta_input.size());
    for (std::size_t i = 0; i < theta.size(); ++i) {
      theta[i] = theta_input[i];
    }

    if(cfg.containsElementNamed("beta")) {
    std::vector<double> beta_input = get_numvec_copy(cfg, "beta");
    beta = CppAD::vector<double>(beta_input.size());
    for (std::size_t i = 0; i < theta.size(); ++i) {
      beta[i] = beta_input[i];
    }
    } else {
      beta = CppAD::vector<double>(0);
    }

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

    if(cfg.containsElementNamed("group")) {
      group_sparsity = cfg["group"];
    } else {
      group_sparsity = Rcpp::List();
    }



  }
};
#endif