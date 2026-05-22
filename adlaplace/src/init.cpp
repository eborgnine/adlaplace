#include <Rcpp.h>

#include "adlaplace/api/density_registry.hpp"
#include "adlaplace/densities/neg_binom.hpp"
#include "adlaplace/densities/random.hpp"

namespace {

using RandomDiagListFn = CppAD::vector<CppAD::AD<double>> (*)(
  const CppAD::vector<CppAD::AD<double>>&,
  const Rcpp::List&,
  const Rcpp::List&);

RandomDiagListFn random_diagonal_list_fn() {
  return static_cast<RandomDiagListFn>(
    static_cast<CppAD::vector<CppAD::AD<double>> (*)(
      const CppAD::vector<CppAD::AD<double>>&,
      const Rcpp::List&,
      const Rcpp::List&)>(&random_diagonal));
}

}  // namespace

void register_adlaplace_default_densities() {
  static bool registered = false;
  if (registered) return;
  registered = true;

  register_log_dens_obs("neg_binom_obs", &neg_binom_obs);
  register_log_dens_single_data("neg_binom_extra", &neg_binom_extra);
  register_log_dens_single_random_diag("random_diagonal", random_diagonal_list_fn());
  register_log_dens_single_random_diag("random", random_diagonal_list_fn());
}
