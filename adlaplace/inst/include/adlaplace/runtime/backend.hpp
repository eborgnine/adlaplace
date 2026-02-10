#ifndef ADLAPLACE_BACKEND_HPP
#define ADLAPLACE_BACKEND_HPP

#include <cppad/cppad.hpp>

struct HessianPack {
  std::vector<int> hessian_p;
  std::vector<int> hessian_i;
  std::vector<int> dim;

  std::vector<int> map_p;
  std::vector<int> map_local;
  std::vector<int> map_global;

};


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

  CPPAD_TESTVECTOR(double) w, wthree, direction_zeros, direction;
  CppAD::sparse_rc<CPPAD_TESTVECTOR(size_t)> unused_pattern;

  CPPAD_TESTVECTOR(double) x;

};


// backend private context
struct BackendContext {
  std::vector<GroupPack> *adPack;
//  std::array<std::array<int,3>,3> sizes; // beta, gamma, theta; begin end N

  size_t Nparams;
  size_t Ngroups;
  size_t Nbeta;
  size_t Ngamma;
  size_t Ntheta;

  HessianPack hessian_inner; // p, local index, global index, dims
  HessianPack hessian_outer; // 

};


#endif
