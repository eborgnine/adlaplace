
#include "adlaplace/adlaplace.hpp"

void gradDense(
	const CppAD::vector<double> &x,
	double &f,
	std::vector<double> &g, 
	std::vector<GroupPack>& adPack
	) {

	const std::size_t n = static_cast<std::size_t>(x.size());

	double fOut = 0.0;
	std::vector<double> gOut(n,0.0);
	CppAD::vector<double> w(1);
	w[0] = 1.0;

	for(size_t D=0;D<adPack.size();D++) {
		CppAD::vector<double> y = adPack[D].fun.Forward(0, x);
		fOut += y[0];
		CppAD::vector<double> gradHere = adPack[D].fun.Reverse(1, w);
		for(size_t Dpar=0;Dpar<n;Dpar++) {
			gOut[Dpar] += gradHere[Dpar];
		}
	}

	f = fOut;
	for(size_t Dpar=0;Dpar<n;Dpar++) {
		g[Dpar] = gOut[Dpar];
	}
}

#include "configObj.hpp"
