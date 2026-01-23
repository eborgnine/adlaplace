#pragma once

#ifndef ADLAPLACE_HPP
#define ADLAPLACE_HPP

// We do NOT support ColPack in this package.#ifndef CPPAD_HAS_COLPACK
#define CPPAD_HAS_COLPACK 0

#include <cppad/cppad.hpp>

// needs cppad but not Rcpp
#include "adlaplace/cppadUtils.hpp"
#include "adlaplace/logexp.hpp"

#include "adlaplace/adpack.hpp"
#include "adlaplace/foromp.hpp"

// needs Eigen, adpack
#include <Eigen/SparseCore>   
#include "adlaplace/functions.hpp"


// for what follows need rcpp

#include <Rcpp.h>
#include "adlaplace/trustOptimUtils.hpp"
// need Eigen 
#include "adlaplace/matrixUtils.hpp"
// needs matrixUtils 
#include "adlaplace/data.hpp"
#include "adlaplace/config.hpp"


// needs adpack
#include "adlaplace/sparsity.hpp"

// needs adpack, omp, config
#include "adlaplace/third.hpp"

// needs config, data, sparsity
#include "adlaplace/adfun.hpp"

// needs adfun
#include "adlaplace/debugging.hpp"


// from trustOptim
#include <CG-sparse.h> 

// needs functions.hpp and trustoptim
#include "adlaplace/innerOpt.hpp"

// needs everything
#include "adlaplace/Rinterfaces_backend.hpp"


#endif
