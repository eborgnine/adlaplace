#ifndef ADLAPLACE_RUNTIME_THREAD_GROUPS_HPP
#define ADLAPLACE_RUNTIME_THREAD_GROUPS_HPP

#include <cstddef>
#include <vector>

#include <Rcpp.h>

#include "adlaplace/api/backend.hpp"
#include "adlaplace/runtime/interfaces_detail.hpp"

inline void adlaplace_require_owner_threads_assigned(const ad_fun& backend) {
  for (std::size_t s = 0; s < backend.fun.size(); ++s) {
    adlaplace_adpack_handle* h = shard_handle(const_cast<ad_fun*>(&backend), s);
    const GroupPack* pack = pack_ctx(h->ctx);
    if (!pack->owner_thread_assigned) {
      Rcpp::stop(
        "OpenMP thread assignment missing; call ad_fun(..., num_threads=) before inner_opt or trace_hinv_t"
      );
    }
  }
}

inline std::vector<std::vector<std::size_t>> thread_groups_from_backend(ad_fun& backend) {
  const std::size_t n = backend.fun.size();
  std::vector<std::size_t> owners(n);
  std::size_t max_t = 0;
  for (std::size_t s = 0; s < n; ++s) {
    const std::size_t t = pack_ctx(backend.fun[s]->ctx)->owner_thread;
    owners[s] = t;
    if (t > max_t) max_t = t;
  }
  std::vector<std::vector<std::size_t>> groups(max_t + 1);
  for (std::size_t s = 0; s < n; ++s) {
    groups[owners[s]].push_back(s);
  }
  return groups;
}

#endif
