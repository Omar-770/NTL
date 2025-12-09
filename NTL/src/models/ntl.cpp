#include "ntl.h"

namespace NTL
{

	double calculate_Z(double Z0, double d, const std::vector<double>& Cn, int M, double z)
	{
		double Z{};
		double temp = 2 * M_PI * z / d;
		

		int N = Cn.size();
		for (int n = 0; n < N - M; n++)
			Z += Cn[n] * std::cos(temp * n);

		int j = N - M - 1;
		for (int m = 1; m < M + 1; m++)
			Z += Cn[j + m] * std::sin(temp * m);

		return Z0 * std::exp(Z);
	}

	double calculate_Z(const NTL& ntl, double z)
	{
		return calculate_Z(ntl.get_Z0(), ntl.get_d(), ntl.get_Cn(), ntl.get_M(), z);
	}

	double calculate_Z(const double* Cn, const size_t& n, const size_t& m,const double& Z0, const double& d, const double& z)
	{
		double Z{};
		double temp = 2 * M_PI * z / d;
		for (size_t i = 0; i < n - m; ++i) {
			Z += Cn[i] * std::cos(temp * i);
		}

		int j = n - m - 1;
		for (size_t i = 1; i < m + 1; ++i) {
			Z += Cn[j + i] * std::cos(temp * i);
		}
		return Z0 * std::exp(Z);
	}

	std::complex<double> calculate_Zin(double Z0, double e_r, double d, const std::vector<double>& Cn, int M, double Zl, double f, int K)
	{
		matrix2x2cd T = calculate_T_matrix(Z0, e_r, d, Cn, M, f, K);

		std::complex<double> A = T(0, 0), B = T(0, 1), C = T(1, 0), D = T(1, 1);
		std::complex<double> Zin = (Zl * A + B) / (Zl * C + D);

		return Zin;
	}

	std::complex<double> calculate_Zin(const NTL& ntl, double Zl, double f, int K)
	{
		return calculate_Zin(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), ntl.get_M(), Zl, f, K);
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

	double calculate_er_eff(double z, double e_r) //function of impedance z
	{
		double wh = calculate_W_H(z, e_r);
		return (e_r + 1) / 2 + (e_r - 1) / 2 * ((std::pow(1 + 12 / wh, -0.5) + ((wh < 1) ? 0.04 * (std::pow(1 - wh, 2)) : 0)));
	}

	matrix2x2cd calculate_T_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, int M, double f, int K)
	{
		double _dz = d / K;		
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

		if(M == 0) // Even optimisation
		{
			matrix2x2cd T_half = matrix2x2cd::Identity();
			int half_K = K / 2;
			for (int i = 0; i < half_K; i++)
			{
				double z = calculate_Z(Z0, d, Cn, M, (double(i) + 0.5) * _dz);
				calculate_Ti(temp_Ti, z, f, _dz, calculate_er_eff(z, e_r));
				T_half *= temp_Ti;
			}

			matrix2x2cd T_half_rev;
			T_half_rev << T_half(1, 1), T_half(0, 1), T_half(1, 0), T_half(0, 0);

			if (K % 2 == 0)
			{
				return T_half * T_half_rev;
			}
			else
			{
				double z_mid = calculate_Z(Z0, d, Cn, M, (double(half_K) + 0.5) * _dz);
				calculate_Ti(temp_Ti, z_mid, f, _dz, calculate_er_eff(z_mid, e_r));
				return T_half * temp_Ti * T_half_rev;
			}
		} 
		else
		{
			matrix2x2cd T = matrix2x2cd::Identity();
			for (int i = 0; i < K; i++)
			{
				double z = calculate_Z(Z0, d, Cn, M, (double(i) + 0.5) * _dz);
				calculate_Ti(temp_Ti, z, f, _dz, calculate_er_eff(z, e_r));
				T *= temp_Ti;
			}

			return T;
		}
	}

	matrix2x2cd calculate_T_matrix(const NTL& ntl, double f, int K)
	{
		return calculate_T_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), ntl.get_M(), f, K);
	}

	std::pair<matrix2x2cd, std::vector<matrix2x2cd>> calculate_T_matrix_with_grad(
		double Z0, double er, double d, const std::vector<double>& Cn, int M, double f, int K)
	{
		double _dz = d / K;
		int N = Cn.size();
		std::vector<matrix2x2cd> L(K), R(K), T(K), dT_dZ(K), dT_dCn(N, matrix2x2cd::Zero());
		std::vector<double> Z(K);

		auto calculate_Ti = [](matrix2x2cd& T, const double& Z, const double& f, const double& d, const double& e_r) {
			const double theta{ 2 * M_PI * f * std::sqrt(e_r) * d / M_C };
			const double cos_theta{ std::cos(theta) };
			const double sin_theta{ std::sin(theta) };

			T(0, 0) = { cos_theta, 0 };
			T(0, 1) = { 0, Z * sin_theta };
			T(1, 0) = { 0, sin_theta / Z };
			T(1, 1) = { cos_theta, 0 };
			};

		matrix2x2cd temp_T;
		double dZ = 1e-6;

		if(M == 0) // Even optimisaton
		{
			int half_K = (K + 1) / 2;
			
			// Calculate
			for (int i = 0; i < half_K; i++)
			{
				Z[i] = calculate_Z(Z0, d, Cn, M, (double(i) + 0.5) * _dz);
				double e_eff = calculate_er_eff(Z[i], er);
				calculate_Ti(T[i], Z[i], f, _dz, e_eff);
				calculate_Ti(temp_T, Z[i] + dZ, f, _dz, e_eff);
				dT_dZ[i] = (temp_T - T[i]) / dZ;
			}

			// Mirror
			for (int i = half_K; i < K; i++)
			{
				T[i] = T[K - 1 - i];
				dT_dZ[i] = dT_dZ[K - 1 - i];
				Z[i] = Z[K - 1 - i];
			}

			// Calculate
			L[0] = matrix2x2cd::Identity();
			for (int i = 1; i < K; i++)
				L[i] = L[i - 1] * T[i - 1];

			// Mirror
			for (int i = 0; i < K; i++)
			{
				R[K - 1 - i] << L[i](1, 1), L[i](0, 1), L[i](1, 0), L[i](0, 0);
			}			
					
		}
		else
		{
			// Calculate
			for (int i = 0; i < K; i++)
			{
				Z[i] = calculate_Z(Z0, d, Cn, M, (double(i) + 0.5) * _dz);
				double e_eff = calculate_er_eff(Z[i], er);
				calculate_Ti(T[i], Z[i], f, _dz, e_eff);
				calculate_Ti(temp_T, Z[i] + dZ, f, _dz, e_eff);
				dT_dZ[i] = (temp_T - T[i]) / dZ;
			}

			// Calculate
			L[0] = matrix2x2cd::Identity();
			R[K - 1] = matrix2x2cd::Identity();
			for (int i = 1; i < K; i++)
			{
				L[i] = L[i - 1] * T[i - 1];
				R[K - 1 - i] = T[K - i] * R[K - i];
			}
			
		}

		// Final 
		for (int k = 0; k < K; k++)
		{
			double zi = (double(k) + 0.5) * _dz;
			for (int i = 0; i < N; i++)
				dT_dCn[i] += L[k] * dT_dZ[k] * R[k] * Z[k] * (i < N - M ? std::cos(2 * M_PI * zi / d * double(i))
																		: std::sin(2 * M_PI * zi / d * double(i - N + M + 1)));
		}

		return { L[K - 1] * T[K - 1], dT_dCn };
	}

	matrix2x2cd calculate_S_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, int M,
		double f, double Zs, double Zl, int K)
	{
		matrix2x2cd T = calculate_T_matrix(Z0, e_r, d, Cn, M, f, K);
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

	matrix2x2cd calculate_S_matrix(const NTL& ntl, double f, double Zs, double Zl, int K)
	{
		return calculate_S_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), ntl.get_M(), f, Zs, Zl, K);
	}

	matrix2x2cd calculate_Y_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, int M,
		double f, int K)
	{
		matrix2x2cd T = calculate_T_matrix(Z0, e_r, d, Cn, M, f, K);
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

	matrix2x2cd calculate_Y_matrix(const NTL& ntl, double f, int K)
	{
		return calculate_Y_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), ntl.get_M(), f, K);
	}

	/// CLASS_NTL_BEGIN

	double NTL::Z(double z) const //function of position
	{

		return calculate_Z(m_Z0, m_d, m_Cn, m_M, z);

	}

	std::complex<double> NTL::Zin(double Zl, double f, int K) const
	{
		return calculate_Zin(*this, Zl, f, K);
	}

	double NTL::W_H(double z) const //function of position z
	{
		double impedance_at_z = Z(z);
		return calculate_W_H(impedance_at_z, m_er);
	}

	double NTL::er_eff(const double& z) const //function of position z
	{
		double impedance_at_z = Z(z);
		return calculate_er_eff(impedance_at_z, m_er);
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

	DATA NTL::data() const
	{
		return DATA{ m_Z0, m_er, m_d, m_Cn };
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