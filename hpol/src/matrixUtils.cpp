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
  const std::vector<std::vector<double>>& randomHessian, 
  const std::vector<double>& qHessian, 
  const Rcpp::List& sparsity, 
  const Config& config, 
  const bool onlyRandom) {

    const size_t Ngroups = randomHessian.size();
      // TO DO, sum sparse matrix
    const Rcpp::List secondAll = config.sparsity["second"];
    const Rcpp::List outAll = onlyRandom?secondAll["random"]:secondAll["full"];
    const Rcpp::IntegerVector iAll = outAll["i"];
    const auto sizeH = iAll.size();
    Rcpp::NumericVector hessianSum(sizeH);


    for(size_t D=0;D < sizeH; D++) {
      hessianSum[D] = 0;
    }

    for(size_t Dgroup = 0;Dgroup<Ngroups;++Dgroup) {


      const Rcpp::List sparsityHere = sparsity[Dgroup];
      const Rcpp::List secondHere = sparsityHere["second"];
      const Rcpp::List targetHere= onlyRandom?secondHere["random"]:secondHere["full"];
      const Rcpp::IntegerVector matchHere = targetHere["match"];

      const std::vector<double> hessianOutHere = randomHessian[Dgroup];
      const size_t Nhere = matchHere.size();
      for(size_t D=0;D < Nhere; D++) {
        hessianSum[matchHere[D]] += hessianOutHere[D];
      }
    } // groupd

      const Rcpp::List secondHere = config.sparsity["Q"];
      const Rcpp::List targetHere= onlyRandom?secondHere["random"]:secondHere["full"];
      const Rcpp::IntegerVector matchHere = targetHere["match"];
      const size_t Nq = matchHere.size();

    for(size_t D=0;D<Nq;++D) {
      hessianSum[matchHere[D]] += qHessian[D];
    }      

    const Rcpp::List outList = onlyRandom?secondAll["random"]:secondAll["full"];
    const Rcpp::IntegerVector outI = outList["i"];
    const Rcpp::IntegerVector outJ = outList["j"];
    const Rcpp::IntegerVector outP = outList["p"];
    const size_t N = outP.size()-1;

    Rcpp::S4 result = make_convert_gCmatrix(hessianSum, outI, outJ, N);

    return(result);
   
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

