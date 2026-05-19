#ifndef HANDLE_CREATE_HPP
#define HANDLE_CREATE_HPP

#include "adlaplace/runtime/backend.hpp"

static void backend_destroy(void* vctx) {
  delete static_cast<AdGroups*>(vctx);
}

#endif
