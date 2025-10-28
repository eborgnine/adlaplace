#pragma once


#include <Rcpp.h>
#include <cppad/cppad.hpp>


#include <type_traits>

#include <utility>      // for std::pair

#include"dataconfig.hpp"
#include"lgamma.hpp"
#include"adfun.hpp"
#include"foromp.hpp"
#include"loglikHelpers.hpp"
#include"matrixUtils.hpp"


#ifndef HPOL_HPP
#define HPOL_HPP


// Constants
inline constexpr double LOGTWOPI      = 1.8378770664093454836;
inline constexpr double HALFLOGTWOPI  = 0.9189385332046727;


#endif // LOGSPACE_HPOL_HPP


