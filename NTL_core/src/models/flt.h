#pragma once

#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <utility>
#include <vector>

#include "models/ntl.h"

#ifdef M_PI
#undef M_PI
#endif 

namespace flt
{
	class flt;

	class flt : public NTL::NTL
	{
	public:
		std::vector<double> _3dB_freqs(double f1, double f2) const;
	};
}