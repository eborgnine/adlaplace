// inst/include/adlaplace/adpack_handle.h
#ifndef ADLAPLACE_ADPACK_HANDLE_H
#define ADLAPLACE_ADPACK_HANDLE_H

#include "adpack_api.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct adlaplace_adpack_handle {
  const adlaplace_adpack_api* api;
  void* ctx;
} adlaplace_adpack_handle;

#ifdef __cplusplus
}
#endif

#endif