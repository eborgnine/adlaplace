#include <Rcpp.h>
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
