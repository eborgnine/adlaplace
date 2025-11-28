
#ifndef DATACONFIG_HPP
#define DATACONFIG_HPP

#include<Rcpp.h>

#include <vector>
#include <cmath>
#include "matrixUtils.hpp"


inline DgCView makeGroups(const Rcpp::List& data) {
  // If a groups matrix is supplied, use it
  if (data.containsElementNamed("groups") && !Rf_isNull(data["groups"])) {
    return DgCView(Rcpp::S4(data["groups"]));
  }

  // Otherwise, build an identity matrix with ncol(XTp) × ncol(XTp)
  const Rcpp::NumericVector y = data["y"];
  const std::size_t Ny = static_cast<std::size_t>(y.size());  // ncol

  return identityMatrix(Ny);
}

struct Data {
  DgCView A;
  DgCView X;
  DgCView groups;

  Rcpp::NumericVector Qdiag, y;
  Rcpp::IntegerVector map;
  size_t Nmap, Nbeta, Ngamma, Ny, Ngroups;

  explicit Data(const Rcpp::List& data)
    : A(         DgCView(Rcpp::S4(data["ATp"])) )
    , X(         DgCView(Rcpp::S4(data["XTp"])) )
    , groups(    makeGroups(data))
    , Qdiag(     data["Qdiag"] )
    , y(         data["y"] )
    , map(       data["map"] )
  {


    Nmap      = static_cast<std::size_t>(map.size());
    Nbeta   = static_cast<std::size_t>(X.nrow());
    Ngamma  = static_cast<std::size_t>(A.nrow());
    Ny    = static_cast<std::size_t>(X.ncol());   // == A.ncol()
    Ngroups = static_cast<std::size_t>(groups.ncol());
    if(Ny != groups.nrow()) {
      Rcpp::Rcout << "Ny " << Ny << " rows of groups " << groups.nrow() << "\n";
    }
    if(Ny != A.ncol()) {
      Rcpp::Rcout << "Ny " << Ny << " columns of A " << A.ncol() << "\n";
    }
    if(Ny != y.size()) {
      Rcpp::Rcout << "lengh y " << y.size() << " columns of X " << Ny << "\n";
    }
  }
};



#endif