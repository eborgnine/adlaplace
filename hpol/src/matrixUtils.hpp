Rcpp::S4 make_dgTMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N, size_t maxEntries);


Rcpp::IntegerVector  compute_p_vector(
  const Rcpp::IntegerVector & j, int ncol);
