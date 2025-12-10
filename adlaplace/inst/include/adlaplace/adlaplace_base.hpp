

#include<Rcpp.h>
#include <cppad/cppad.hpp>
#include "adlaplace/foromp.hpp"
#include <RcppEigen.h>

// needs cppad
#include "adlaplace/cppadUtils.hpp"
#include "adlaplace/logexp.hpp"

// need Rcpp
#include "adlaplace/matrixUtils.hpp"

// needs matrixUtils 
#include "adlaplace/data.hpp"
#include "adlaplace/config.hpp"

