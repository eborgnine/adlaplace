#ifndef LOGSPACE_ADD_HPP
#define LOGSPACE_ADD_HPP



class atomic_logspace_add : public CppAD::atomic_base<double> {
public:
    atomic_logspace_add(const std::string& name);
    
    
  // Forward
    virtual bool forward(
        size_t p,
        size_t q,
        const CppAD::vector<bool>& vx,
        CppAD::vector<bool>& vy,
        const CppAD::vector<double>& tx,
        CppAD::vector<double>& ty
    );

    // Reverse
    virtual bool reverse(
        size_t q,
        const CppAD::vector<double>& tx,
        const CppAD::vector<double>& ty,
        CppAD::vector<double>& px,
        const CppAD::vector<double>& py
    );
};



// Global atomic instance
extern atomic_logspace_add logspace_add_atomic;

// ======= Template for ALL autodiff types =======
template <class Type>
inline Type logspace_add_n(const CppAD::vector<Type>& x) {
    CppAD::vector<Type> yout(1);
    logspace_add_atomic(x, yout);
    return yout[0];
}


template <>
double logspace_add_n(const CppAD::vector<double>& x);


#endif // LOGSPACEADD_HPP