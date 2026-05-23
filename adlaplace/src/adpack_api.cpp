#include "adlaplace/eval/fgh.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"
#include "adlaplace/api/register.hpp"

extern "C" const adlaplace_adpack_api adlaplace_AD_API = {
  ADLAPLACE_ADPACK_API_VERSION,
  1,
  &eval_f,
  &eval_grad,
  &eval_hess,
  &get_sizes,
  &get_sparse_pattern,
  &get_hessian,
  &eval_trace_hinv_t,
  &backend_destroy,
  NULL
};
