#ifndef LOGDENSEXTRAEMPTY_HPP
#define LOGDENSEXTRAEMPTY_HPP

// dummy funciton if there is no logDensExtra
template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	) {
	CppAD::vector<Type> result(1);
	result[0] = 0.0;
	return(result);
}

// declare
template CppAD::vector<double>
logDensExtra<double>(
    const CppAD::vector<double>&,
    const Data&, const Config&);

template CppAD::vector<CppAD::AD<double>>
logDensExtra<CppAD::AD<double>>(
    const CppAD::vector<CppAD::AD<double>>&,
    const Data&, const Config&);


#endif