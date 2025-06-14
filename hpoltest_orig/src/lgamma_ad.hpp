#ifndef LGAMMA_AD_HPP
#define LGAMMA_AD_HPP


//#include <cppad/cppad.hpp>


class atomic_lgamma_ad : public CppAD::atomic_base<double> {
public:
  atomic_lgamma_ad(const std::string& name);

virtual bool forward(
    size_t p,
    size_t q,
    const CppAD::vector<bool>& vx,
    CppAD::vector<bool>& vy,
    const CppAD::vector<double>& tx,
    CppAD::vector<double>& ty
) override;

virtual bool reverse(
    size_t q,
    const CppAD::vector<double>& tx,
    const CppAD::vector<double>& ty,
    CppAD::vector<double>& px,
    const CppAD::vector<double>& py
) override;
};

// Declaration of the global atomic function instance
extern atomic_lgamma_ad lgamma_ad_atomic;

inline double lgamma_ad(double x) {
    return std::lgamma(x);
}

template<class T>
inline CppAD::AD<T> lgamma_ad(CppAD::AD<T> x) {
    return CppAD::AD<T>(lgamma_ad(Value(x)));
}
#endif
