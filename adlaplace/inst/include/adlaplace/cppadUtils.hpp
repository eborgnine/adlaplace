#ifndef CPPADUTILS_HPP
#define CPPADUTILS_HPP

#include <cppad/cppad.hpp>


inline CppAD::vector<CppAD::AD<double>> slice(
  const CppAD::vector<CppAD::AD<double>>& x,
  size_t start, size_t end)
{
    CppAD::vector<CppAD::AD<double>> out(end - start);
    for (size_t i = start, j = 0; i < end; ++i, ++j)
        out[j] = x[i];
    return out;
}

inline double to_double(double x) {
    return x;
}

inline double to_double(const CppAD::AD<double>& x) {
    return CppAD::Value(x);
}

#endif