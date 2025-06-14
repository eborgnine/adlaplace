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


// Declaration of the global atomic function instance
extern atomic_logspace_add logspace_add_atomic;

// logspace_add for doubles (safe, just use std::log(std::exp(x) + std::exp(y)))
inline double logspace_add(double x, double y) {
    if (x == -INFINITY) return y;
    if (y == -INFINITY) return x;
    if (x > y)
        return x + std::log1p(std::exp(y - x));
    else
        return y + std::log1p(std::exp(x - y));
}

// logspace_add for AD<T>
template<class T>
inline CppAD::AD<T> logspace_add(CppAD::AD<T> x, CppAD::AD<T> y) {
    // fallback if you want an overloaded version
    return CppAD::AD<T>(logspace_add(Value(x), Value(y)));
}

// You probably want a vectorized version as well:
CppAD::AD<double> logspace_add_n(const CppAD::vector<CppAD::AD<double>>& x);

template<class Type>
Type logspace_add_n(const CppAD::vector<Type>& x);

#endif // LOGSPACEADD_HPP