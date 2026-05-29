#include "adlaplace/runtime/ad_fun_clone.hpp"

#include "adlaplace/api/register.hpp"

GroupPack clone_group_pack(const GroupPack& src) {
  GroupPack dst;
  dst.fun = src.fun;
  dst.work_grad = src.work_grad;
  dst.work_hess = src.work_hess;
  dst.work_inner_grad = src.work_inner_grad;
  dst.work_inner_hess = src.work_inner_hess;
  dst.pattern_grad = src.pattern_grad;
  dst.pattern_grad_inner = src.pattern_grad_inner;
  dst.pattern_hessian = src.pattern_hessian;
  dst.pattern_hessian_inner = src.pattern_hessian_inner;
  dst.x = src.x;
  dst.w = src.w;
  dst.wthree = src.wthree;
  dst.direction_zeros = src.direction_zeros;
  dst.direction = src.direction;
  dst.unused_pattern = src.unused_pattern;
  dst.shard_index = src.shard_index;
  dst.owner_thread = src.owner_thread;
  dst.n_beta = src.n_beta;
  dst.n_theta = src.n_theta;
  return dst;
}

ad_fun* clone_ad_fun(const ad_fun* src) {
  if (!src) {
    Rcpp::stop("clone_ad_fun: source handle is NULL");
  }
  if (src->fun.empty()) {
    Rcpp::stop("clone_ad_fun: source has no AD shards");
  }

  std::size_t n_beta = 0;
  std::size_t n_theta = 0;
  std::vector<GroupPack> packs;
  packs.reserve(src->fun.size());

  for (adlaplace_adpack_handle* h : src->fun) {
    if (!h || !h->ctx) {
      continue;
    }
    const GroupPack* pack = static_cast<const GroupPack*>(h->ctx);
    if (packs.empty()) {
      n_beta = pack->n_beta;
      n_theta = pack->n_theta;
    }
    packs.push_back(clone_group_pack(*pack));
  }

  if (packs.empty()) {
    Rcpp::stop("clone_ad_fun: source has no valid shards");
  }

  return packs_to_ad_fun(std::move(packs), n_beta, n_theta);
}
