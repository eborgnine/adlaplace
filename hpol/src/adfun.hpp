#pragma once
#include"dataconfig.hpp"



struct GroupPack {
  CppAD::ADFun<double>                    fun;       // taped function for the group
  CppAD::sparse_hessian_work              work;      // reusable work cache
  CppAD::vector< std::set<size_t> >       pattern;   // Hessian sparsity pattern (size = Nparams)
  std::array< std::vector<size_t>, 3 >    outRowCol; // [0]=rows, [1]=cols for subset extraction
  std::vector<size_t> outP;
  std::array< std::vector<size_t>, 3 >    nsRowCol; 
  std::vector<size_t> nsP;
};

struct AdpackHandle {
  std::vector<GroupPack>* ptr = nullptr;   // pointer to existing or new object
  bool created_here = false;             // whether we must delete it later

  void cleanup() {
    if (created_here && ptr) { delete ptr; ptr = nullptr; }
  }
};


// Forward declaration of helper
CPPAD_TESTVECTOR(std::set<size_t>)
build_pattern_from_R(const Rcpp::IntegerVector& row0,
                     const Rcpp::IntegerVector& col0,
                     size_t n);

AdpackHandle getAdpackFromR(SEXP adfun,
                            const std::vector<double>& parametersC,
                            const Data& dataC,
                            const Config& configC);


CppAD::ADFun<double> adFunGroup(
  const std::vector<double> & parameters,  
  const Data& data, 
  const Config& config,
  const Rcpp::IntegerVector& strataI,
  const size_t start,
  const size_t end
  );

// Declared interface
std::vector<GroupPack> getAdFun(const std::vector<double>& parameters,
                                const Data& data,
                                const Config& config);

CppAD::ADFun<double> adFunQ(
  const std::vector<double> & parameters,  
  const Data& data,
  const Config& config);


GroupPack getAdFunQ(
  const std::vector<double>& parameters,
               const Data&                data,
               const Config&              config);