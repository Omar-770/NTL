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

	double calculate_Z(const NTL& ntl, double z)
	{
		return calculate_Z(ntl.get_Z0(), ntl.get_d(), ntl.get_Cn(), z);
	}

	double calculate_Z(const double* Cn, const size_t& n, const double& Z0, const double& d, const double& z)
	{
		double Z{};
		double temp = 2 * M_PI * z / d;
		for (size_t i = 0; i < n; ++i) {
			Z += Cn[i] * std::cos(temp * i);
		}
		return Z0 * std::exp(Z);
	}

	std::complex<double> calculate_Zin(double Z0, double e_r, double d, const std::vector<double>& Cn, double Zl, double f, int K)
	{
		matrix2x2cd T = calculate_T_matrix(Z0, e_r, d, Cn, f, K);

		std::complex<double> A = T(0, 0), B = T(0, 1), C = T(1, 0), D = T(1, 1);
		std::complex<double> Zin = (Zl * A + B) / (Zl * C + D);
		
		return Zin;
	}

	std::complex<double> calculate_Zin(const NTL& ntl, double Zl, double f, int K)
	{
		return calculate_Zin(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), Zl, f, K);
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

	matrix2x2cd calculate_T_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, int K)
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
			calculate_Ti(temp_Ti, z, f, _dz, calculate_er_eff(z, e_r));
			T *= temp_Ti;
		}

		return T;
	}

	matrix2x2cd calculate_T_matrix(const NTL& ntl, double f, int K)
	{
		return calculate_T_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), f, K);
	}

	std::pair<matrix2x2cd, std::vector<matrix2x2cd>> calculate_T_matrix_with_grad(
		double Z0, double er, double d, const std::vector<double>& Cn, double f, int K)
	{
		double _dz = d / K;
		int N = Cn.size();
		auto calculate_Ti = [](matrix2x2cd& T, const double& Z, const double& f, const double& d, const double& e_r) {
			const double theta{ 2 * M_PI * f * std::sqrt(e_r) * d / M_C };
			const double cos_theta{ std::cos(theta) };
			const double sin_theta{ std::sin(theta) };

			T(0, 0) = { cos_theta, 0 };
			T(0, 1) = { 0, Z * sin_theta };
			T(1, 0) = { 0, sin_theta / Z };
			T(1, 1) = { cos_theta, 0 };
			};

		std::vector<matrix2x2cd> L(K), R(K), T(K), dT_dCn(N, matrix2x2cd::Zero());

		for (int i = 0; i < K; i++)
		{
			double z = calculate_Z(Z0, d, Cn, (double(i) + 0.5) * _dz);
			calculate_Ti(T[i], z, f, _dz, er);
		}

		L[0] = matrix2x2cd::Identity();
		R[K - 1] = matrix2x2cd::Identity();

		for (int i = 1; i < K; i++)
		{
			L[i] = L[i - 1] * T[i - 1];
			R[K - 1 - i] = T[K - i] * R[K - i];
		}

		matrix2x2cd temp_T;
		double dZ = 1e-6;

		for (int i = 0; i < N; i++)
			for (int k = 0; k < K; k++)
			{
				double zi = (double(k) + 0.5) * _dz;
				double z = calculate_Z(Z0, d, Cn, zi);
				calculate_Ti(temp_T, z + dZ, f, _dz, er);
				dT_dCn[i] += L[k] * (temp_T - T[k]) / dZ * R[k] * z * std::cos(2 * M_PI * zi / d * double(i));
			}	


		return { L[K - 1] * T[K - 1], dT_dCn };
	}

	
	std::pair<matrix2x2cd, std::vector<matrix2x2cd>> GMN_calculate_T_matrix_with_grad(
		double Z0, double er, double d, const std::vector<double>& Cn, double f, int K)
	{
		double _dz = d / K;
		int N = Cn.size();

		std::vector<matrix2x2cd> Ti(K);
		std::vector<double> Zi(K);
		std::vector<matrix2x2cd> dTi_dZ(K);

		// Helper lambda for local T calculation
		auto calc_Ti_local = [&](double Z, double len, double eps) {
			matrix2x2cd T;
			double e_eff = calculate_er_eff(Z, eps);
			double theta = 2 * M_PI * f * std::sqrt(e_eff) * len / M_C;
			double c = std::cos(theta);
			double s = std::sin(theta);
			T << std::complex<double>(c, 0), std::complex<double>(0, Z * s),
				std::complex<double>(0, s / Z), std::complex<double>(c, 0);
			return T;
			};

		// 1. Pre-calculate Section Matrices and Local Derivatives
		// (Parallelizable loop)
		//#pragma omp parallel for
		for (int i = 0; i < K; i++)
		{
			double z_pos = (double(i) + 0.5) * _dz;
			Zi[i] = calculate_Z(Z0, d, Cn, z_pos);

			Ti[i] = calc_Ti_local(Zi[i], _dz, er);

			// Local finite difference for dTi/dZ (robust)
			double h = Zi[i] * 1e-6;
			matrix2x2cd Ti_plus = calc_Ti_local(Zi[i] + h, _dz, er);
			dTi_dZ[i] = (Ti_plus - Ti[i]) / h;
		}

		// 2. Forward Sweep (Prefix Products)
		std::vector<matrix2x2cd> L(K);
		L[0] = matrix2x2cd::Identity();
		for (int i = 1; i < K; i++)
			L[i] = L[i - 1] * Ti[i - 1];

		// 3. Backward Sweep (Suffix Products)
		std::vector<matrix2x2cd> R(K);
		R[K - 1] = matrix2x2cd::Identity();
		for (int i = K - 2; i >= 0; i--)
			R[i] = Ti[i + 1] * R[i + 1];

		// 4. Calculate Total T
		matrix2x2cd T_total = L[K - 1] * Ti[K - 1];

		// 5. Assemble Gradients
		std::vector<matrix2x2cd> grads(N, matrix2x2cd::Zero());

		// d(Z)/dCn term helper
		auto dZ_dCn = [&](int i, int n) {
			double z_pos = (double(i) + 0.5) * _dz;
			double term = 2 * M_PI * z_pos / d;
			// Derivative of exp( sum(Cn * cos) ) is Z * cos
			return Zi[i] * std::cos(term * n);
			};

		// Combine: dT/dCn = Sum_i ( L_i * dTi/dZ * R_i * dZ_i/dCn )
		//#pragma omp parallel for
		for (int n = 0; n < N; n++)
		{
			for (int i = 0; i < K; i++)
			{
				// M_i = L[i] * dTi_dZ[i] * R[i]
				// This matrix is "Sensitivity of Total T to Impedance Z_i"
				matrix2x2cd M_i = L[i] * dTi_dZ[i] * R[i];
				grads[n] += M_i * dZ_dCn(i, n);
			}
		}

		return { T_total, grads };
	}

	matrix2x2cd calculate_S_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, double f, double Zs, double Zl, int K)
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

	matrix2x2cd calculate_S_matrix(const NTL& ntl, double f, double Zs, double Zl, int K)
	{
		return calculate_S_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), f, Zs, Zl, K);
	}

	matrix2x2cd calculate_Y_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, int K)
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

	matrix2x2cd calculate_Y_matrix(const NTL& ntl, double f, int K)
	{
		return calculate_Y_matrix(ntl.get_Z0(), ntl.get_er(), ntl.get_d(), ntl.get_Cn(), f, K);
	}
	
	/// CLASS_NTL_BEGIN

	double NTL::Z(double z) const //function of position
	{

		return calculate_Z(m_Z0, m_d, m_Cn, z);

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