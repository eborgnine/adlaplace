
#ifndef ADPACK_HPP
#define ADPACK_HPP

#include <cppad/cppad.hpp>

struct GroupPack {
  CppAD::ADFun<double>              fun;       // taped function for the group
  CppAD::sparse_jac_work            work_grad;
  CppAD::sparse_hes_work            work_hess;      // reusable work cache
  CppAD::sparse_jac_work            work_inner_grad;
  CppAD::sparse_hes_work            work_inner_hess;      // reusable work cache
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_grad_inner;
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian;  // note these are upper triangle only
  CppAD::sparse_rcv<CPPAD_TESTVECTOR(size_t), CPPAD_TESTVECTOR(double)> pattern_hessian_inner;

  CPPAD_TESTVECTOR(double) w;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

};


#endif