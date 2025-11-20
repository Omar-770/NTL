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

namespace NTL
{
	class WPD;
	struct WPD_DATA;

	using matrix2x2cd = Eigen::Matrix<std::complex<double>, 2, 2>;
	using matrix3x3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
	using matrix4x4cd = Eigen::Matrix<std::complex<double>, 4, 4>;
	using matrix5x5cd = Eigen::Matrix<std::complex<double>, 5, 5>;
}