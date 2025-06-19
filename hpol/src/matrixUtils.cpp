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




Rcpp::S4 make_dgTMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N, size_t maxEntries,
  bool isLower)
{

  int Ni = N;
  size_t maxEntries2 = maxEntries == 0? 1: maxEntries-1;

  Rcpp::IntegerVector dims = Rcpp::IntegerVector::create(Ni, Ni);
  Rcpp::IntegerVector iR, jR;
  Rcpp::NumericVector xR;


  if(!isLower) {
// Convert std::vector<size_t> to IntegerVector (R requires int32 indices)
      xR = x[Rcpp::Range(0, maxEntries2)];
      iR = i[Rcpp::Range(0, maxEntries2)];
      jR = j[Rcpp::Range(0, maxEntries2)];
  } else {
        // Need to symmetrize: count entries
    int n_in = x.size();
    int n_out = 0;
    for(int k=0; k<n_in; ++k) {
      n_out++; // original entry
      if(i[k] != j[k]) n_out++; // add symmetric if off-diagonal
    }

    xR = Rcpp::NumericVector(n_out);
    iR = Rcpp::IntegerVector(n_out);
    jR = Rcpp::IntegerVector(n_out);

    for (int Dlower = 0,Dfull=0; Dlower < n_in; ++Dlower,++Dfull) {
        xR[Dfull] = x[Dlower];
        iR[Dfull] = i[Dlower];
        jR[Dfull] = j[Dlower];
        if(i[Dlower] != j[Dlower]) { //fill in upper
            Dfull++;
            xR[Dfull] = x[Dlower];
            iR[Dfull] = j[Dlower];
            jR[Dfull] = i[Dlower];
        }   
    }
}


Rcpp::S4 mat("dgTMatrix");
mat.slot("i") = iR;
mat.slot("j") = jR;
mat.slot("x") = xR;
mat.slot("Dim") = dims;
//  mat.slot("uplo") =  Rcpp::CharacterVector uplo= Rcpp::wrap('L');
return mat;



}
