// macros to declare log density functions, for double and CppAD::AD<double>
// used in objectiveFunction.cpp

#ifndef INSTANTIATE_LOGDENS
#define INSTANTIATE_LOGDENS

#define INSTANTIATE_LOGDENSOBS(T) \
template CppAD::vector<CppAD::AD<double>> logDensObs<T>( \
  const CppAD::vector<CppAD::AD<double>>&, \
  const CppAD::vector<T>&, \
  const CppAD::vector<T>&, \
  const Data&, const Config&, const size_t);

#define INSTANTIATE_LOGDENSEXTRA(T) \
template CppAD::vector<T> logDensExtra<T>( \
  const CppAD::vector<T>&, \
  const Data&, const Config&);

#define INSTANTIATE_LOGDENSRANDOM(T) \
template CppAD::vector<CppAD::AD<double>> logDensRandom<T>( \
	const CppAD::vector<CppAD::AD<double>>&, \
	const CppAD::vector<T> &, \
	const Data& data, const Config& config);

INSTANTIATE_LOGDENSOBS(double)
INSTANTIATE_LOGDENSOBS(CppAD::AD<double>)

INSTANTIATE_LOGDENSEXTRA(double)
INSTANTIATE_LOGDENSEXTRA(CppAD::AD<double>)

INSTANTIATE_LOGDENSRANDOM(double)
INSTANTIATE_LOGDENSRANDOM(CppAD::AD<double>)

#undef INSTANTIATE_LOGDENSOBS
#undef INSTANTIATE_LOGDENSEXTRA
#undef INSTANTIATE_LOGDENSRANDOM

#endif