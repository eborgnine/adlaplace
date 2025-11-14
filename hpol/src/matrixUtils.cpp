#include "hpol.hpp"
#include <vector>


Rcpp::IntegerVector compute_p_vector(
    const Rcpp::IntegerVector& j, 
    int ncol) {
    // Step 1: Count nonzeros per column
    std::vector<int> counts(ncol, 0);
    for (size_t k = 0; k < j.size(); ++k) {
        // j[k] is zero-based column index
        counts[j[k]]++;
    }

    // Step 2: Cumulative sum (CSC pointer)   
    Rcpp::IntegerVector p(ncol + 1);
    p[0] = 0;
    for (int c = 0; c < ncol; ++c) {
        p[c + 1] = p[c] + counts[c];
    }
//    return p;std::vector<int> p(ncol + 1, 0);
//    std::partial_sum(counts.begin(), counts.end(), p.begin() + 1);

    return p;
}



Rcpp::S4 make_TMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N)
{

  int Ni = N;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);

  Rcpp::LogicalVector toKeep = j >= 0;
  Rcpp::IntegerVector iR=i[toKeep], jR=j[toKeep];
  Rcpp::NumericVector xR=x[toKeep];


  Rcpp::S4 mat("dsTMatrix");
  mat.slot("i") = iR;
  mat.slot("j") = jR;
  mat.slot("x") = xR;
  mat.slot("Dim") = dims;
  mat.slot("uplo") =  Rcpp::wrap('L');
  return mat;



}

Rcpp::S4 make_CMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& p)
{

  int Ni = p.size()-1;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);

  Rcpp::S4 mat("dsCMatrix");
  mat.slot("i") = i;
  mat.slot("p") = p;
  mat.slot("x") = x;
  mat.slot("Dim") = dims;
  mat.slot("uplo") =  Rcpp::wrap('U');
  return mat;



}

Rcpp::S4 make_gCMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& p)
{

  int Ni = p.size()-1;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);

  Rcpp::S4 mat("dgCMatrix");
  mat.slot("i") = i;
  mat.slot("p") = p;
  mat.slot("x") = x;
  mat.slot("Dim") = dims;
  return mat;

}


Rcpp::S4 make_convert_gCmatrix(const Rcpp::NumericVector& x,
                              const Rcpp::IntegerVector& i,
                              const Rcpp::IntegerVector& j,
                              const size_t N)
{
  const int nnz = i.size();

  std::vector<int> Iout;
  std::vector<int> Jout;
  std::vector<double> Xout;
  Iout.reserve(2 * nnz);
  Jout.reserve(2 * nnz);
  Xout.reserve(2 * nnz);
  for(int k=0;k<nnz;k++) {
    const int row = i[k];
    const int col = j[k];
    const double val = x[k]; 

      Iout.push_back(row);
    Jout.push_back(col);
    Xout.push_back(val);

    // add mirrored lower triangle if off-diagonal
    if (row != col) {
      Iout.push_back(col);
      Jout.push_back(row);
      Xout.push_back(val);
    }
  }
    const int nnzOut = Jout.size();
    std::vector<size_t> idx(nnzOut);
    std::iota(idx.begin(), idx.end(), 0);  

    std::stable_sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) { return Jout[a] < Jout[b]; });

  Rcpp::IntegerVector i1(nnzOut), j1(nnzOut);
  Rcpp::NumericVector x1(nnzOut);    
    for(int D=0;D<nnzOut;D++) {
      const int newD = idx[D];
      i1[D] = Iout[newD];
      j1[D] = Jout[newD];
      x1[D] = Xout[newD];
    }
  Rcpp::IntegerVector p1 = compute_p_vector(j1, N);
  Rcpp::S4 result = make_gCMatrix(x1, i1, p1);
  return(result);
}

Rcpp::S4 assembleHessian(
  std::vector<GroupPack>& adpack, 
  GroupPack& qPack, 
  const Rcpp::List& sparsity, 
  const Config& config, 
  const bool onlyRandom) {

    const size_t Ngroups = adpack.size();

    const Rcpp::List secondAll = config.sparsity["second"];
    const Rcpp::List outAll = onlyRandom?secondAll["random"]:secondAll["full"];
    const Rcpp::IntegerVector iAll = outAll["i"];
    const Rcpp::IntegerVector jAll = outAll["j"];
    const Rcpp::IntegerVector pAll = outAll["p"];

    const int N = static_cast<int>(pAll.size() - 1L);
    const auto sizeH = iAll.size();
    Rcpp::NumericVector hessianSum(sizeH, 0.0);


    //for(size_t D=0;D < sizeH; D++) {hessianSum[D] = 0;}
      if (config.verbose ) Rcpp::Rcout << "groups.. ";

    for(size_t Dgroup = 0;Dgroup<Ngroups;++Dgroup) {

      const std::vector<size_t>& matchHere = adpack[Dgroup].outRowCol[2];
      const size_t Nhere = matchHere.size();
      const CppAD::vector<double>& hessianOutHere = adpack[Dgroup].out_hess.val();

      for(size_t D=0;D < Nhere; D++) {
        const size_t indexHere = matchHere[D];
        hessianSum[indexHere] += hessianOutHere[D];
      }

    } // groupd

      if (config.verbose ) Rcpp::Rcout << "done groups\n";

// Q
      const std::vector<size_t>& matchHere = qPack.outRowCol[2];
      const size_t Nhere = matchHere.size();
      const CppAD::vector<double>& hessianOutHere = qPack.out_hess.val();

      for(size_t D=0;D < Nhere; D++) {
        const size_t indexHere = matchHere[D];
        hessianSum[indexHere] += hessianOutHere[D];
      }
      if (config.verbose ) Rcpp::Rcout << "convert to gCmatrix\n";

    Rcpp::S4 result = make_convert_gCmatrix(hessianSum, iAll, jAll, N);

    return(result);
   
}


