#ifndef DATACONFIG_HPP
#define DATACONFIG_HPP

#include "adlaplace/dgcview.hpp"
#include "adlaplace/cppadUtils.hpp"

struct Data {
  DgCView A;
  DgCView X;

  Rcpp::NumericVector Qdiag, y;
  DgCView map;
  DgCView elgm_matrix;

  size_t Nmap, Nbeta, Ngamma, Ny;

  explicit Data(const Rcpp::List& data)
    : A(         DgCView(Rcpp::S4(data["ATp"])) )
    , X(         DgCView(Rcpp::S4(data["XTp"])) )
    , Qdiag(     data["Qdiag"] )
    , y(         data["y"] )
    , map(       DgCView(Rcpp::S4(data["map"])) )
  {

    if(data.containsElementNamed("elgm_matrix")) {
      elgm_matrix = DgCView(Rcpp::S4(data["elgm_matrix"]));
    } 
    Nmap = map.ncol();
    Nbeta   = static_cast<std::size_t>(X.nrow());
    Ngamma  = static_cast<std::size_t>(A.nrow());
    Ny    = static_cast<std::size_t>(X.ncol());   // == A.ncol()
    if(Ny != A.ncol()) {
      Rcpp::Rcout << "Ny " << Ny << " columns of A " << A.ncol() << "\n";
    }
    if(Ny != y.size()) {
      Rcpp::Rcout << "lengh y " << y.size() << " columns of X " << Ny << "\n";
    }
  }
};

// ----- config bundle -----
struct Config {
  bool verbose;
  bool transform_theta;
  int num_threads=1;

  CppAD::vector<double> beta;
  CppAD::vector<double> start_gamma;
  CppAD::vector<double> theta;    

  DgCView groups;

  std::vector<std::vector<int>> hessianIP_outer = std::vector<std::vector<int>>(2);
  std::vector<std::vector<int>> hessianIPLower_outer = std::vector<std::vector<int>>(3);

  std::vector<std::vector<int>> hessianIP_inner = std::vector<std::vector<int>>(2);
  std::vector<std::vector<int>> hessianIPLower_inner = std::vector<std::vector<int>>(3);

  Rcpp::List group_sparsity_inner;
  Rcpp::List group_sparsity_outer;


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

    if(cfg.containsElementNamed("groups")) {
      groups = DgCView(Rcpp::S4(cfg["groups"]));
    } 

    if(cfg.containsElementNamed("hessian_inner")) {
      const Rcpp::S4 sm = cfg["hessian_inner"];
      const Rcpp::IntegerVector i = sm.slot("i");
      const Rcpp::IntegerVector p = sm.slot("p");
      hessianIP_inner[0] = std::vector<int>(i.begin(), i.end());
      hessianIP_inner[1] = std::vector<int>(p.begin(), p.end());
    }
    if(cfg.containsElementNamed("hessian_outer")) {
      const Rcpp::S4 sm = cfg["hessian_outer"];
      const Rcpp::IntegerVector i = sm.slot("i");
      const Rcpp::IntegerVector p = sm.slot("p");
      hessianIP_outer[0] = std::vector<int>(i.begin(), i.end());
      hessianIP_outer[1] = std::vector<int>(p.begin(), p.end());
    }

    if(cfg.containsElementNamed("hessianL_inner")) {
      const Rcpp::S4 sm = cfg["hessianL_inner"];
      const Rcpp::IntegerVector i = sm.slot("i");
      const Rcpp::IntegerVector p = sm.slot("p");
      const Rcpp::IntegerVector x = sm.slot("x");
      hessianIPLower_inner[0] = std::vector<int>(i.begin(), i.end());
      hessianIPLower_inner[1] = std::vector<int>(p.begin(), p.end());
      hessianIPLower_inner[2] = std::vector<int>(x.begin(), x.end());
    }
    if(cfg.containsElementNamed("hessianL_outer")) {
      const Rcpp::S4 sm = cfg["hessianL_outer"];
      const Rcpp::IntegerVector i = sm.slot("i");
      const Rcpp::IntegerVector p = sm.slot("p");
      const Rcpp::IntegerVector x = sm.slot("x");
      hessianIPLower_outer[0] = std::vector<int>(i.begin(), i.end());
      hessianIPLower_outer[1] = std::vector<int>(p.begin(), p.end());
      hessianIPLower_outer[2] = std::vector<int>(x.begin(), x.end());
    }

    if(cfg.containsElementNamed("group_outer")) {
      group_sparsity_outer = cfg["group_outer"];
    } else {
      group_sparsity_outer = Rcpp::List();
    }
    if(cfg.containsElementNamed("group_inner")) {
      group_sparsity_inner = cfg["group_inner"];
    } else {
      group_sparsity_inner = Rcpp::List();
    }

  }
};

#endif
