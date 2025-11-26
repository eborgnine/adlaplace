#ifndef ADFUN_HPP
#define ADFUN_HPP
#include"dataconfig.hpp"

struct GroupPack {
  CppAD::ADFun<double>                    fun;       // taped function for the group
  CppAD::sparse_hes_work              work_hess;      // reusable work cache
  CppAD::sparse_rc< CppAD::vector<size_t> >  pattern_hess;   // Hessian sparsity pattern (size = Nparams)
  CppAD::sparse_rcv< CppAD::vector<size_t>, CppAD::vector<double> > out_hess; // values
  CppAD::sparse_jac_work work_grad;
  CppAD::sparse_rc< CppAD::vector<size_t> > pattern_grad;
  CppAD::sparse_rcv< CppAD::vector<size_t>, CppAD::vector<double> > out_grad;
  std::array< std::vector<size_t>, 3 >    nsRowCol; 
  std::vector<size_t> nsP;
  std::array< std::vector<size_t>, 3 >    outRowCol; // [0]=rows, [1]=cols for subset extraction
  std::vector<size_t> outP;
};

struct AdpackHandle {
  std::vector<GroupPack>* ptr = nullptr;   // pointer to existing or new object
  bool created_here = false;             // whether we must delete it later

  void cleanup() {
    if (created_here && ptr) { delete ptr; ptr = nullptr; }
  }
};






#endif
