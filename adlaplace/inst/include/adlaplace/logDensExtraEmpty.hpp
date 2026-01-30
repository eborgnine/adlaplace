#ifndef LOGDENSEXTRAEMPTY_HPP
#define LOGDENSEXTRAEMPTY_HPP

// dummy funciton if there is no logDensExtra

CppAD::vector<CppAD::AD<double>> logDensExtra(
	const CppAD::vector<CppAD::AD<double>> &theta,
	const Data& data,
	const Config& config
	) {
	CppAD::vector<CppAD::AD<double>> result(1);
	result[0] = 0.0;
	return(result);
}

#endif
