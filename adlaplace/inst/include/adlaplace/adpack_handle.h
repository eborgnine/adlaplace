#ifndef ADLAPLACE_ADPACK_HANDLE_H
#define ADLAPLACE_ADPACK_HANDLE_H

#include<stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Versioning so you can evolve the API safely
#define ADLAPLACE_ADPACK_API_VERSION 1

typedef struct adlaplace_adpack_api {
  int api_version;     // must be ADLAPLACE_ADPACK_API_VERSION
  int thread_safe;     // 1 if safe under OpenMP (no R calls)

  // Evaluate f, grad, hess for "inner" problem (gamma), at x (length Ngamma).
  // Return 0 on success; nonzero error code otherwise.
  int (*f)(void* ctx, const int *i, const double* x, double* out_f);
  int (*f_grad)(void* ctx, const int *i, const double* x, const bool *inner, 
    double* out_f, double* out_grad);
  int (*f_grad_hess)(void* ctx, const int *i, const double* x, const bool *inner,
                 double* out_f,
                     double* out_grad, double* out_hess);


int (*get_sizes)(void* ctx, size_t* Nparams, size_t* Ngroups,
                 size_t* Nbeta, size_t* Ngamma, size_t* Ntheta);


int (*get_hessian)(void* ctx, const bool *inner,
                            const int** p, size_t* p_len,
                            const int** i, size_t* i_len);


  // Optional destructor for ctx
  void (*destroy)(void* ctx);

  // Optional: return last error message (thread-local if parallel)
  const char* (*last_error)(void* ctx);

} adlaplace_adpack_api;

typedef struct adlaplace_adpack_handle {
  const adlaplace_adpack_api* api;
  void* ctx;

} adlaplace_adpack_handle;




#ifdef __cplusplus
}
#endif

#endif