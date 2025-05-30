#ifndef LOGSPACE_HPOL_HPP
#define LOGSPACE_HPOL_HPP



#include <Rcpp.h>
#include <cppad/cppad.hpp>

using CppAD::AD;
using CppAD::vector;

// Constants
#define LOGTWOPI 1.8378770664093454835606594728112352797227949472755668
#define HALFLOGTWOPI 0.918938533204672741780329736405617639861397473637783

#include"lgamma.hpp"
#include"logspaceadd.hpp"


#endif