#ifndef CPPADUTILS_HPP
#define CPPADUTILS_HPP

#include <cppad/cppad.hpp>

static const std::string JAC_COLOR = "cppad";  
static const std::string HESS_COLOR = "cppad.symmetric";

struct GroupPack {
  CppAD::ADFun<double>              fun;       // taped function for the group
  CppAD::sparse_jac_work            work_grad;
  CppAD::sparse_hes_work            work_hess;      // reusable work cache
};


struct AdpackHandle {
  std::vector<GroupPack>* ptr = nullptr;   // pointer to existing or new object
  bool created_here = false;             // whether we must delete it later

  void cleanup() {
    if (created_here && ptr) { delete ptr; ptr = nullptr; }
  }
};


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
