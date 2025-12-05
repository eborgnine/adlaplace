
#ifndef DATACONFIG_HPP
#define DATACONFIG_HPP


//#include <vector>
//#include <cmath>


struct Data {
  DgCView A;
  DgCView X;

  Rcpp::NumericVector Qdiag, y;
  Rcpp::IntegerVector map;

  DgCView elgm_matrix;

  size_t Nmap, Nbeta, Ngamma, Ny;

  explicit Data(const Rcpp::List& data)
    : A(         DgCView(Rcpp::S4(data["ATp"])) )
    , X(         DgCView(Rcpp::S4(data["XTp"])) )
    , Qdiag(     data["Qdiag"] )
    , y(         data["y"] )
    , map(       data["map"] )
  {

    if(data.containsElementNamed("elgm_matrix")) {
      elgm_matrix = DgCView(Rcpp::S4(data["elgm_matrix"]));
    } 
    Nmap      = static_cast<std::size_t>(map.size());
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



#endif