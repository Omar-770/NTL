#include "ntl.h"

namespace NTL
{

	double calculate_Z(double Z0, double d, const std::vector<double>& Cn, double z)
	{
		double Z{};
		double temp = 2 * M_PI * z / d;

		for (int n = 0; n < Cn.size(); n++)
			Z += Cn[n] * std::cos(temp * n);

		return Z0 * std::exp(Z);
	}

	double calculate_Z(const NTL& ntl, const double z)
	{
		return calculate_Z(ntl.get_Z0(), ntl.get_d(), ntl.get_Cn(), z);
	}

	double calculate_Z(const double* Cn, size_t n, const double& Z0, const double& d, const double& z)
	{
		double Z{};
		double temp = 2 * M_PI * z / d;
		for (size_t i = 0; i < n; ++i) {
			Z += Cn[i] * std::cos(temp * i);
		}
		return Z0 * std::exp(Z);
	}


	double calculate_W_H(double z, double e_r) // z is IMPEDANCE
	{
		double e_eff_cross = (e_r + 1) / 2.0 + (e_r - 1) / 2.0 * std::pow(1.0 + 12.0 / 2.0, -0.5);
		double Z_cross = (60.0 / std::sqrt(e_eff_cross)) * std::log(8.0 / 2.0 + 2.0 / 4.0);

		if (z >= Z_cross)
		{
			double A = (z / 60) * std::sqrt((e_r + 1) / 2) + ((e_r - 1) / (e_r + 1)) * (0.23 + (0.11 / e_r));
			return (8 * std::exp(A)) / (std::exp(2 * A) - 2);
		}
		else
		{
			double B = (377 * M_PI) / (2 * z * std::sqrt(e_r));
			return (2 / M_PI) * (B - 1 - std::log(std::abs(2 * B - 1)) + ((e_r - 1) / (2 * e_r)) * (std::log(std::abs(B - 1)) + 0.39 - (0.61 / e_r)));
		}
	}

	double calculate_e_eff(double z, double e_r) //function of impedance z
	{
		double wh = calculate_W_H(z, e_r);
		return (e_r + 1) / 2 + (e_r - 1) / 2 * ((std::pow(1 + 12 / wh, -0.5) + ((wh < 1) ? 0.04 * (std::pow(1 - wh, 2)) : 0)));
	}

	matrix2x2cd calculate_T_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, double K)
	{
		double _dz = d / K;
		matrix2x2cd T = matrix2x2cd::Identity();
		matrix2x2cd temp_Ti;

		auto calculate_Ti = [](matrix2x2cd& T, const double& Z, const double& f, const double& d, const double& e_r) {
			const double theta{ 2 * M_PI * f * std::sqrt(e_r) * d / M_C };
			const double cos_theta{ std::cos(theta) };
			const double sin_theta{ std::sin(theta) };

			T(0, 0) = { cos_theta, 0 };
			T(0, 1) = { 0, Z * sin_theta };
			T(1, 0) = { 0, sin_theta / Z };
			T(1, 1) = { cos_theta, 0 };
			};

		for (int i = 0; i < K; i++)
		{
			double z = calculate_Z(Z0, d, Cn, (double(i) + 0.5) * _dz);
			calculate_Ti(temp_Ti, z, f, _dz, calculate_e_eff(z, e_r));
			T *= temp_Ti;
		}

		return T;
	}

	matrix2x2cd calculate_T_matrix(const NTL& ntl, double f, double K)
	{
		return calculate_T_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), f, K);
	}

	matrix2x2cd calculate_S_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, double f, double Zs, double Zl, double K)
	{
		matrix2x2cd T = calculate_T_matrix(Z0, e_r, d, Cn, f, K);
		std::complex<double> A = T(0, 0);
		std::complex<double> B = T(0, 1);
		std::complex<double> C = T(1, 0);
		std::complex<double> D = T(1, 1);

		matrix2x2cd S;

		std::complex<double> denominator = A * Zl + B + C * Zs * Zl + D * Zs;

		if (std::abs(denominator) < 1e-12)
		{
			S(0, 0) = 0; S(0, 1) = 0;
			S(1, 0) = 0; S(1, 1) = 0;
			return S;
		}

		S(0, 0) = (A * Zl + B - C * Zs * Zl - D * Zs) / denominator;
		S(0, 1) = (2.0 * (A * D - B * C) * std::sqrt(Zs * Zl)) / denominator;
		S(1, 0) = (2.0 * std::sqrt(Zs * Zl)) / denominator;
		S(1, 1) = (-A * Zl + B - C * Zs * Zl + D * Zs) / denominator;

		return S;
	}

	matrix2x2cd calculate_S_matrix(const NTL& ntl, double f, double Zs, double Zl, double K)
	{
		return calculate_S_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), f, Zs, Zl, K);
	}

	matrix2x2cd calculate_Y_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, double K)
	{
		matrix2x2cd T = calculate_T_matrix(Z0, e_r, d, Cn, f, K);
		std::complex<double> A = T(0, 0);
		std::complex<double> B = T(0, 1);
		std::complex<double> C = T(1, 0);
		std::complex<double> D = T(1, 1);

		matrix2x2cd Y;
		std::complex<double> denominator = B;
		if (std::abs(denominator) < 1e-12)
		{
			Y(0, 0) = 0; Y(0, 1) = 0;
			Y(1, 0) = 0; Y(1, 1) = 0;
			return Y;
		}
		Y(0, 0) = D / denominator;
		Y(0, 1) = -(A * D - B * C) / denominator;
		Y(1, 0) = -1.0 / denominator;
		Y(1, 1) = A / denominator;
		return Y;
	}

	matrix2x2cd calculate_Y_matrix(const NTL& ntl, double f, double K)
	{
		return calculate_Y_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), f, K);
	}
	
	/// CLASS_NTL_BEGIN

	double NTL::Z(double z) const //function of position
	{

		return calculate_Z(m_Z0, m_d, m_Cn, z);

	}

	double NTL::W_H(double z) const //function of position z
	{
		double impedance_at_z = Z(z);
		return calculate_W_H(impedance_at_z, m_er);
	}

	double NTL::e_eff(const double& z) const //function of position z
	{
		double impedance_at_z = Z(z);
		return calculate_e_eff(impedance_at_z, m_er);
	}

	std::vector<std::pair<double, double>> NTL::get_Z_vec(double step_size) const
	{
		std::vector<std::pair<double, double>> Z_vec;
		Z_vec.reserve(m_d / step_size + 1);
		for (double z = 0; z <= m_d; z += step_size)
		{
			Z_vec.emplace_back(z, Z(z));
		}
		return Z_vec;
	}

	std::vector<std::pair<double, double>> NTL::get_w_h_vec(double step_size) const
	{
		std::vector<std::pair<double, double>> w_h_vec;
		w_h_vec.reserve(m_d / step_size + 1);

		for (double z = 0; z <= m_d; z += step_size)
		{
			w_h_vec.emplace_back(z, W_H(z));
		}

		return w_h_vec;
	}

	NTL_DATA NTL::data() const
	{
		return NTL_DATA{ m_Z0, m_er, m_d, m_Cn };
	}

	matrix2x2cd NTL::T_matrix(double f, int K) const
	{
		return calculate_T_matrix(*this, f, K);
	}

	matrix2x2cd NTL::S_matrix(double f, double Zs, double Zl, int K) const
	{
		return calculate_S_matrix(*this, f, Zs, Zl, K);
	}

	matrix2x2cd NTL::Y_matrix(double f, int K) const
	{
		return calculate_Y_matrix(*this, f, K);
	}


	std::complex<double> NTL::S11(double f, double Zs, double Zl, int K) const
	{
		return calculate_S_matrix(*this, f, Zs, Zl, K)(0, 0);
	}

	/// CLASS_NTL_END
}