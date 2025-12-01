#ifndef ADFUN_HPP
#define ADFUN_HPP

#include "adlaplace/adpack.hpp"
#include "adlaplace/cppadUtils.hpp"
#include "adlaplace/sparsity.hpp"

template <class Type>
CppAD::vector<CppAD::AD<double>> logDensObs(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &beta,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config,
	const size_t Dgroup
	);

template <class Type>
CppAD::vector<Type> logDensExtra(
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	);

template <class Type>
CppAD::vector<CppAD::AD<double>> logDensRandom(
	const CppAD::vector<CppAD::AD<double>>& gamma,
	const CppAD::vector<Type> &theta,
	const Data& data,
	const Config& config
	);

std::vector<GroupPack> getAdFunOuter(
	const Data& data,
	const Config& config) {

	std::vector<GroupPack> result(data.Ngroups+2L);

	if(config.verbose) {
		Rcpp::Rcout << " groups " << data.Ngroups << " ";
	}

		const size_t Ngamma = config.start_gamma.size();
		const size_t Nbeta = config.beta.size();
		const size_t Ntheta = config.theta.size();
		const size_t Nparams = Nbeta + Ngamma + Ntheta;
		const size_t nbSizeIndex = Nparams-1;

# pragma omp parallel
	{
		CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
		for(size_t D=0;D<Nbeta;D++) {
			ad_params[D] = config.beta[D];
		}
		for(size_t Dgamma=0;Dgamma<Ngamma;Dgamma++) {
			ad_params[Dgamma+Nbeta] = config.start_gamma[Dgamma];
		}
		for(size_t Dtheta=0;Dtheta<Ntheta;Dtheta++) {
				ad_params[Dtheta+Nbeta+Ngamma] = config.theta[Dtheta];
		} 


# pragma omp for nowait
		for(size_t D=0;D<data.Ngroups;D++) {

			CppAD::Independent(ad_params);
			auto resultHere = logDensObs(
				slice(ad_params, Nbeta, Nbeta + Ngamma), // gamma
				slice(ad_params, 0, Nbeta), // beta
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config, D);
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);
		}



# pragma omp single nowait
		{
			CppAD::Independent(ad_params);

			auto resultHere = logDensRandom(
				slice(ad_params, Nbeta, Nbeta + Ngamma), // gamma
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[data.Ngroups].fun = std::move(fun);

		}

		# pragma omp single nowait
		{
			CppAD::Independent(ad_params);

			auto resultHere = logDensExtra(
				slice(ad_params, Nbeta + Ngamma, Nparams), // theta
				data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[data.Ngroups+1].fun = std::move(fun);

		}

	} // parallel

	if(config.group_sparsity.size()) {
		if(config.verbose) {
			Rcpp::Rcout << "add sparse\n";
		}

		addSparsityToAdFun(result, config.group_sparsity);

	}

	return(result);
}


std::vector<GroupPack> getAdFunInner(
	const Data& data,
	const Config& config) {

	auto parameters = config.start_gamma;


	std::vector<GroupPack> result(data.Ngroups+1L);

	if(config.verbose) {
		Rcpp::Rcout << " adfun Nparams " << parameters.size() << " groups " << data.Ngroups << " ";
	}

# pragma omp parallel
	{
		const size_t Nparams = parameters.size();
		CppAD::vector<CppAD::AD<double>> ad_params(Nparams);
		for(size_t D=0;D<Nparams;D++) {
			ad_params[D] = parameters[D];
		}


# pragma omp for nowait
		for(size_t D=0;D<data.Ngroups;D++) {

			CppAD::Independent(ad_params);   
			auto resultHere = logDensObs(
				ad_params, 
				config.beta, 
				config.theta,
				data, config, D);
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[D].fun = std::move(fun);
		}



# pragma omp single nowait
		{
			CppAD::Independent(ad_params);
			auto resultHere = logDensRandom(ad_params, 
				config.theta, data, config);   
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[data.Ngroups].fun = std::move(fun);

		}

		# pragma omp single nowait
		{
			// adds the extra part to the likelihood
			// won't contribute to derivatives
			CppAD::Independent(ad_params);
			CppAD::vector<double> resultExtra = logDensExtra(
				config.theta, 
				data, config);   
			CppAD::vector<CppAD::AD<double>> resultHere(1);
			resultHere[0] = resultExtra[0];
			CppAD::ADFun<double> fun(ad_params, resultHere);
			result[data.Ngroups+1].fun = std::move(fun);

		}

	} // parallel

	if(config.group_sparsity.size()) {
		if(config.verbose) {
			Rcpp::Rcout << "add sparse\n";
		}

		addSparsityToAdFun(result, config.group_sparsity);

	}

	return(result);
}


#endif
