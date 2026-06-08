#ifndef ADLAPLACE_EXTENSION_BACKEND_ADPACK_API_HPP
#define ADLAPLACE_EXTENSION_BACKEND_ADPACK_API_HPP

#include "adlaplace/api/adpack_handle.h"
#include "adlaplace/eval/assign_memory.hpp"
#include "adlaplace/eval/fgh.hpp"
#include "adlaplace/eval/trace_hinv_t.hpp"

void backend_destroy(void* vctx);

#define ADLAPLACE_DEFINE_BACKEND_ADPACK_API(api_symbol) \
  extern "C" const adlaplace_adpack_api api_symbol = { \
    ADLAPLACE_ADPACK_API_VERSION, \
    1, \
    &eval_f, \
    &eval_grad, \
    &eval_hess, \
    &get_sizes, \
    &get_sparse_pattern, \
    &get_hessian, \
    &eval_assign_memory, \
    &eval_trace_hinv_t, \
    &backend_destroy, \
    NULL \
  };

#endif
