#ifndef ADLAPLACE_AD_FUN_CLONE_HPP
#define ADLAPLACE_AD_FUN_CLONE_HPP

#include "adlaplace/api/backend.hpp"

GroupPack clone_group_pack(const GroupPack& src);

ad_fun* clone_ad_fun(const ad_fun* src);

#endif
