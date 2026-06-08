#ifndef ADLAPLACE_EVAL_ASSIGN_MEMORY_IMPL_HPP
#define ADLAPLACE_EVAL_ASSIGN_MEMORY_IMPL_HPP

#include "adlaplace/eval/assign_memory.hpp"

#include <algorithm>
#include <cppad/cppad.hpp>
#include <cstddef>

#include "adlaplace/api/backend.hpp"

int eval_assign_memory(void* vctx) {
  if (vctx == nullptr) return 1;

  GroupPack& gp = *pack_ctx(vctx);
  const std::size_t n_params = gp.fun.Domain();
  if (n_params == 0) return 2;

  gp.direction.clear();
  gp.direction_zeros.clear();
  gp.wthree.clear();

  gp.direction.resize(n_params);
  gp.direction_zeros.resize(n_params);
  std::fill(gp.direction_zeros.begin(), gp.direction_zeros.end(), 0.0);
  gp.wthree.resize(3);
  gp.wthree[0] = 0.0;
  gp.wthree[1] = 0.0;
  gp.wthree[2] = 1.0;

  gp.fun.capacity_order(3);

  return 0;
}

#endif
