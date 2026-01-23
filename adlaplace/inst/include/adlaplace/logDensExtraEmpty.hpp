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

#endif
