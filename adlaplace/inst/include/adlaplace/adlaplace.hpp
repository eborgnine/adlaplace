#ifndef ADLAPLACE_HPP
#define ADLAPLACE_HPP


#include "adlaplace/adlaplace_base.hpp"

//  only needs Rcpp, RcppEigen, cppad
#include "adlaplace/trustOptimUtils.hpp"

// needs data, config
#include "adlaplace/adpack.hpp"

// needs adpack
#include "adlaplace/sparsity.hpp"
#include "adlaplace/functions.hpp"
#include "adlaplace/adfun.hpp"

// needs sparsity
#include "adlaplace/third.hpp"

// needs functions.hpp, stuff in trustOptim package
#include "adlaplace/innerOpt.hpp"

// needs everything
#include "adlaplace/Rinterfaces.hpp"

// not done lgamma 


#endif
