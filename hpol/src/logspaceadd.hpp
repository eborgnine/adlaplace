#ifndef LOGSPACE_ADD_HPP
#define LOGSPACE_ADD_HPP

#include <cppad/cppad.hpp>


class atomic_logspace_add : public CppAD::atomic_four<double> {
public:
    size_t n;
    atomic_logspace_add(const std::string& name, size_t n_);
    
    bool for_type(
        size_t call_id,
        const CppAD::vector<CppAD::ad_type_enum>& type_x,
        CppAD::vector<CppAD::ad_type_enum>& type_y
    ) override;
    
    bool forward(
        size_t call_id,
        const CppAD::vector<bool>& select_y,
        size_t order_low,
        size_t order_up,
        const CppAD::vector<double>& tx,
        CppAD::vector<double>& ty
    ) override;
    
    bool reverse(
        size_t call_id,
        const CppAD::vector<bool>& select_x,
        size_t order_up,
        const CppAD::vector<double>& tx,
        const CppAD::vector<double>& ty,
        CppAD::vector<double>& px,
        const CppAD::vector<double>& py
    ) override;
};



template<class Type>
Type logspace_add_ad(const CppAD::vector<Type>& x, atomic_logspace_add& atomic) {
    CppAD::vector<Type> y(1);
    atomic(x, y); // call the atomic instance
    return y[0];
}



#endif