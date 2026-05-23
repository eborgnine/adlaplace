#include <Rcpp.h>

#include "adlaplace/api/density_registry.hpp"
#include "adlaplace/densities/neg_binom.hpp"
#include "adlaplace/densities/random.hpp"

void register_adlaplace_default_densities() {
  static bool registered = false;
  if (registered) return;
  registered = true;

  register_log_dens_obs("neg_binom_obs", &neg_binom_obs);
  register_log_dens_single_data("neg_binom_extra", &neg_binom_extra);
  register_log_dens_single_random_diag("random_diagonal", &random_diagonal);
  register_log_dens_single_random_diag("random", &random_diagonal);
}
