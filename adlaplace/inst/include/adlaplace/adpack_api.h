// inst/include/adlaplace/adpack_api.h
#ifndef ADLAPLACE_ADPACK_API_H
#define ADLAPLACE_ADPACK_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Versioning so you can evolve the API safely
#define ADLAPLACE_ADPACK_API_VERSION 1

typedef struct adlaplace_adpack_api {
  int api_version;     // must be ADLAPLACE_ADPACK_API_VERSION
  int thread_safe;     // 1 if safe under OpenMP (no R calls)
  size_t Ngamma;

  // Evaluate f, grad, hess for "inner" problem (gamma), at x (length Ngamma).
  // Return 0 on success; nonzero error code otherwise.
  int (*f)(void* ctx, const int *i, const double* x, double* out_f);
  int (*f_grad)(void* ctx, const int *i, const double* x, double* out_f, double* out_grad);
  int (*f_grad_hess)(void* ctx, const int *i, const double* x, double* out_f,
                     double* out_grad, double* out_hess);

  // Optional destructor for ctx
  void (*destroy)(void* ctx);

  // Optional: return last error message (thread-local if parallel)
  const char* (*last_error)(void* ctx);

} adlaplace_adpack_api;

#ifdef __cplusplus
} // extern "C"
#endif

#endif