#ifndef HANDLE_CREATE_HPP
#define HANDLE_CREATE_HPP

#include "adlaplace/runtime/backend.hpp"

static void backend_destroy(void* vctx) {
  delete static_cast<GroupPack*>(vctx);
}

static void ad_groups_destroy(ad_groups* groups) {
  if (!groups) return;
  for (adlaplace_adpack_handle* h : groups->fun) {
    if (!h) continue;
    if (h->api && h->api->destroy && h->ctx) {
      h->api->destroy(h->ctx);
    }
    delete h;
  }
  delete groups;
}

#endif
