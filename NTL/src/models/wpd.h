#pragma once

#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <utility>
#include <vector>
#include <omp.h>

#include "models/ntl.h"

#ifdef M_PI
#undef M_PI
#endif 

namespace WPD
{
	class WPD;
	struct WPD_DATA;

	using matrix2x2cd = Eigen::Matrix<std::complex<double>, 2, 2>;
	using matrix3x3cd = Eigen::Matrix<std::complex<double>, 3, 3>;
	using matrix4x4cd = Eigen::Matrix<std::complex<double>, 4, 4>;
	using matrix5x5cd = Eigen::Matrix<std::complex<double>, 5, 5>;

	inline constexpr double M_PI = 3.14159265358979311599796346854;
	inline constexpr double M_C = 299792458;


	matrix3x3cd calculate_Y_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2,
		double d3, const std::vector<double>& Cn3, double R, double f, int K = 50);
	matrix3x3cd calculate_Y_matrix(const WPD& wpd, double f, int K = 50);

	matrix3x3cd calculate_S_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2,
		double d3, const std::vector<double>& Cn3, double R, double f, std::array<std::complex<double>, 3> Zl, int K = 50);
	matrix3x3cd calculate_S_matrix(const WPD& wpd, double f, std::array<std::complex<double>, 3> Zl, int K = 50);


	struct WPD_DATA
	{
		NTL::NTL ntl2, ntl3;
		double R;
	};

	class WPD
	{
	public:
		WPD() : m_ntl2(), m_ntl3(), m_R(0)
		{

		}

		WPD(const NTL::NTL& ntl2, const NTL::NTL& ntl3, double R)
			: m_ntl2(ntl2), m_ntl3(ntl3), m_R(R), m_Z0(ntl2.get_Z0()), m_er(ntl2.get_er())
		{
			if (ntl2.get_Z0() != ntl3.get_Z0() || ntl2.get_er() != ntl3.get_er())
				throw(std::logic_error("Attempting to initialise a WPD using two different substrates"));
		}

		WPD(const WPD_DATA& data)
			:m_ntl2(data.ntl2), m_ntl3(data.ntl3), m_R(data.R), m_Z0(data.ntl2.get_Z0()), m_er(data.ntl2.get_er())
		{
			if (data.ntl2.get_Z0() != data.ntl3.get_Z0() || data.ntl2.get_er() != data.ntl3.get_er())
				throw(std::logic_error("Attempting to initialise a WPD using two different substrates"));
		}

		matrix3x3cd Y_matrix(double f, int K = 50) const;
		matrix3x3cd S_matrix(double f, std::array<std::complex<double>, 3> Zl, int K = 50) const;

	private:
		NTL::NTL m_ntl2, m_ntl3;
		double m_R;
		double m_Z0;
		double m_er;

	public:
		//setters and getters
		void set_R(const double& R) { m_R = R; };
		void set_arms(const NTL::NTL& ntl2, const NTL::NTL& ntl3) {if (ntl2.get_Z0() != ntl3.get_Z0() || ntl2.get_er() != ntl3.get_er())
				throw(std::logic_error("Attempting to set a WPD using two different substrates")); m_ntl2 = ntl2; m_ntl3 = ntl3; }

		NTL::NTL get_ntl2() const { return m_ntl2; }
		NTL::NTL get_ntl3() const { return m_ntl3; }
		double get_R() const { return m_R; }
		double get_Z0() const { return m_Z0; }
		double get_er() const { return m_er; }
	};
}