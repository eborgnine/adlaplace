#ifndef ADLAPLACE_EXTENSION_BACKEND_EVAL_AMALGAM_HPP
#define ADLAPLACE_EXTENSION_BACKEND_EVAL_AMALGAM_HPP

// Include in exactly one .cpp per shared library (adlaplace.so or a backend .so).
// Defines eval_f, eval_grad, eval_hess, get_sizes, eval_assign_memory,
// eval_trace_hinv_t, and related symbols referenced by backend_adpack_api.hpp.

#include "adlaplace/eval/sizes_impl.hpp"
#include "adlaplace/eval/fgh_impl.hpp"
#include "adlaplace/eval/assign_memory_impl.hpp"
#include "adlaplace/eval/trace_hinv_t_eval_impl.hpp"

#endif
