#pragma once

#ifndef LOGEXP_HPP
#define LOGEXP_HPP


// declare primary template
template<class T>
inline T exp_any(T x);

// specializations (MUST be in a header)
template<>
inline double exp_any<double>(double x) {
  return std::exp(x);
}

template<>
inline CppAD::AD<double> exp_any< CppAD::AD<double> >(CppAD::AD<double> x) {
  return CppAD::exp(x);
}

// declare primary template
template<class T>
inline T log_any(T x);

// specializations (MUST be in a header)
template<>
inline double log_any<double>(double x) {
  return std::log(x);
}

template<>
inline CppAD::AD<double> log_any< CppAD::AD<double> >(CppAD::AD<double> x) {
  return CppAD::log(x);
}

#endif
