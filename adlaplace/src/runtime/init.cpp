#include <Rcpp.h>

#include "adlaplace/runtime/thread_affinity_debug.hpp"
#define ADLAPLACE_MATH_LGAMMA_DEFINE
#include "adlaplace/math/lgamma.hpp"
#define ADLAPLACE_MATH_LOG_ERFC_DEFINE
#include "adlaplace/math/log_erfc.hpp"

void adlaplace_init_atomics() {
  static bool initialized = false;
  if (initialized) return;
  initialized = true;

  adlaplace_debug_print_load_banner();

  adlaplace_init_lgamma_atomic();
  adlaplace_init_log_erfc_atomic();
}
