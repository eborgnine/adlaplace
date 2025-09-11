Rcpp::S4 make_TMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N);


Rcpp::IntegerVector  compute_p_vector(
  const Rcpp::IntegerVector & j, int ncol);
