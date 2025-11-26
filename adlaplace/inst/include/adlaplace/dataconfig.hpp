
#ifndef DATACONFIG_HPP
#define DATACONFIG_HPP

#include<Rcpp.h>

#include <vector>
#include <cmath>
#include "matrixUtils.hpp"


struct Data {
  DgCView A;
  DgCView X;

  Rcpp::NumericVector Qdiag, y;
  Rcpp::IntegerVector map;
  size_t Nmap, Nbeta, Ngamma, Ny;

  explicit Data(const Rcpp::List& data)
    : A(         DgCView(Rcpp::S4(data["ATp"])) )
    , X(         DgCView(Rcpp::S4(data["XTp"])) )
    , Qdiag(     data["Qdiag"] )
    , y(         data["y"] )
    , map(       data["map"] )
  {
    Nmap      = static_cast<std::size_t>(map.size());
    Nbeta   = static_cast<std::size_t>(X.nrow());
    Ngamma  = static_cast<std::size_t>(A.nrow());
    Ny    = static_cast<std::size_t>(X.ncol());   // == A.ncol()
  }
};



#endif