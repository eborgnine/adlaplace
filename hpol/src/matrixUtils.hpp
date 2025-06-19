Rcpp::S4 make_dgTMatrix(
  const Rcpp::NumericVector& x,
  const Rcpp::IntegerVector& i,
  const Rcpp::IntegerVector& j,
  size_t N, size_t maxEntries);


std::vector<int> compute_p_vector(
  const std::vector<int>& j, int ncol);
