#include <Rcpp.h>
#include <vector>
#include <numeric>  // for std::partial_sum

std::vector<int> compute_p_vector(const std::vector<int>& j, int ncol) {
    // Step 1: Count nonzeros per column
    std::vector<int> counts(ncol, 0);
    for (size_t k = 0; k < j.size(); ++k) {
        // j[k] is zero-based column index
        counts[j[k]]++;
    }

    // Step 2: Cumulative sum (CSC pointer)
    std::vector<int> p(ncol + 1, 0);
    std::partial_sum(counts.begin(), counts.end(), p.begin() + 1);

    // p[0] = 0, p[ncol] = total nnz
    return p;
}




Rcpp::S4 make_dgTMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N, size_t maxEntries)
{

  int Ni = N;
  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);
  size_t maxEntries2 = maxEntries == 0? 1: maxEntries-1;
  
    // Convert std::vector<size_t> to IntegerVector (R requires int32 indices)
  Rcpp::NumericVector xR = x[Rcpp::Range(0, maxEntries2)];
  Rcpp::IntegerVector iR = i[Rcpp::Range(0, maxEntries2)];
  Rcpp::IntegerVector jR = j[Rcpp::Range(0, maxEntries2)];


  Rcpp::S4 mat("dgTMatrix");
  mat.slot("i") = iR;
  mat.slot("j") = jR;
  mat.slot("x") = xR;
  mat.slot("Dim") = dims;
//  mat.slot("uplo") =  Rcpp::CharacterVector uplo= Rcpp::wrap('L');
  return mat;
}
