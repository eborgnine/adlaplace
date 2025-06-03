
#include <cppad/cppad.hpp>
#include <set>
#include <vector>
#include <iostream>

int main() {
    using CppAD::AD;
    size_t Nparams = 3;
    CppAD::vector<AD<double>> ad_params(Nparams);
    for (size_t i = 0; i < Nparams; ++i)
        ad_params[i] = i + 1.0;

    CppAD::Independent(ad_params);

    CppAD::vector<AD<double>> y(1);
    y[0] = ad_params[0] * ad_params[0];

    CppAD::ADFun<double> f(ad_params, y);

    CppAD::vector< std::set<size_t> > select_range(1);
    select_range[0].insert(0);
    bool transpose = false;

    auto hes_pattern = f.RevSparseHes(Nparams, select_range, transpose);
    std::cout << "hes_pattern.size() = " << hes_pattern.size() << std::endl;
    for (size_t i = 0; i < hes_pattern.size(); ++i) {
        std::cout << "i=" << i << ":";
        for (auto j : hes_pattern[i])
            std::cout << " " << j;
        std::cout << std::endl;
    }
}

