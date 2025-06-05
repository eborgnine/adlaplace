#ifndef LOGSPACE_ADD_HPP
#define LOGSPACE_ADD_HPP

#include <cppad/cppad.hpp>


class atomic_logspace_add : public CppAD::atomic_four<double> {
public:
    atomic_logspace_add(const std::string& name);
    
private:
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

// Declaration of the global atomic function instance
extern atomic_logspace_add logspace_add_instance;

// Template function declaration
template<class Type>
Type logspace_add_ad(Type x, Type y);

// Explicit instantiation for common types
extern template CppAD::AD<double> logspace_add_ad<CppAD::AD<double>>(CppAD::AD<double>, CppAD::AD<double>);



#endif