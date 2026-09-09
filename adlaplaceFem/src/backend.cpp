// Backend boilerplate for adlaplaceFem (FEM Matérn random densities).
// Taping and evaluation must live in this DSO (macOS CppAD requirement).
// register_impl.hpp provides packs_to_ad_fun / make_ad_pack_ptr (no adlaplace.so link).

#include <Rcpp.h>

#include "adlaplace/extension.hpp"

#include "adlaplace/eval_impl.hpp"
#include "adlaplace/register_impl.hpp"

// Must follow eval_impl.hpp: FemSsqShard reuses pack_sparsity_sizes/get_pattern.
#include "fem_ssq_shard.hpp"

#include <memory>
#include <utility>
#include <vector>

ADLAPLACE_DEFINE_BACKEND(adlaplace_fem_make_shard)

// packs_to_ad_fun's ShardFactory signature is ad_shard*(AdTape&&), which
// cannot carry the analytic payload, so the ssq shard gets its own builder
// rather than a global registry the factory would have to pop from.
ADLAPLACE_EXTENSION_EXPORT ad_pack *
fem_ssq_packs_to_ad_fun(std::vector<AdTape> &&packs, std::size_t n_beta,
                        std::size_t n_theta,
                        std::shared_ptr<const femssq::Data> data) {
  auto *groups = new ad_pack();
  groups->abi_version = ADLAPLACE_ABI_VERSION;
  groups->fun.reserve(packs.size());
  for (std::size_t g = 0; g < packs.size(); ++g) {
    packs[g].shard_index = g;
    packs[g].n_beta = n_beta;
    packs[g].n_theta = n_theta;
    groups->fun.push_back(new FemSsqShard(std::move(packs[g]), data));
  }
  return groups;
}
