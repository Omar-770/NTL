#pragma once

#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <utility>
#include <vector>

#ifdef M_PI
#undef M_PI
#endif 

namespace NTL
{
	class NTL;
	struct NTL_DATA;

	using matrix2x2cd = Eigen::Matrix<std::complex<double>, 2, 2>;

	inline constexpr double M_PI = 3.14159265358979311599796346854;
	inline constexpr double M_C = 299792458;


	double calculate_Z(double Z0, double d, const std::vector<double>& Cn, double z);
	double calculate_Z(const NTL& ntl, const double z);
	double calculate_Z(const double* Cn, size_t n, const double& Z0, const double& d, const double& z); //for nlopt

	double calculate_W_H(double z, double e_r);
	double calculate_e_eff(double z, double e_r);

	matrix2x2cd calculate_T_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, double K = 50);
	matrix2x2cd calculate_T_matrix(const NTL& ntl, double f, double K = 50);

	matrix2x2cd calculate_S_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, double Zs, double Zl, double K = 50);
	matrix2x2cd calculate_S_matrix(const NTL& ntl, double f, double Zs, double Zl, double K = 50);

	matrix2x2cd calculate_Y_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn,
		double f, double K = 50);
	matrix2x2cd calculate_Y_matrix(const NTL& ntl, double f, double K = 50);


	struct NTL_DATA
	{
		double Z0;
		double er;
		double d;

		std::vector<double> Cn;
	};

	class NTL
	{
	public:
		NTL() : m_Z0(0), m_er(0), m_d(0), m_Cn({ 0 }) {

		}

		NTL(double Z0, double er, double d)
			: m_Z0(Z0), m_er(er), m_d(d) {
			m_Cn = std::vector<double>{ 0.0 };
		};

		NTL(double Z0, double er, double d, const std::vector<double>& Cn)
			: m_Z0(Z0), m_er(er), m_d(d), m_Cn(Cn) {
		};


		NTL(const NTL_DATA& data)
			: m_Z0(data.Z0), m_er(data.er), m_d(data.d), m_Cn(data.Cn) {
		};

		//impedance profiles
		double Z(double z) const;
		double W_H(double z) const;
		double e_eff(const double& z) const;

		//data extraction
		std::vector<std::pair<double, double>> get_Z_vec(double step_size = 0.1) const;
		std::vector<std::pair<double, double>> get_w_h_vec(double step_size = 0.1) const;

		NTL_DATA data() const;

		//microwave matrices
		matrix2x2cd T_matrix(double f, int K = 50) const;
		matrix2x2cd S_matrix(double f, double Zs, double Zl, int K = 50) const;
		matrix2x2cd Y_matrix(double f, int K = 50) const;

		//reflection coefficient
		std::complex<double> S11(double f, double Zs, double Zl, int K = 50) const;

	private:

		//Physical parameters
		double m_Z0; // Characteristic impedance
		double m_er;	// Relative permittivity
		double m_d;  // Length of the line

		std::vector<double> m_Cn; // Coefficients vector

	public:

		//setters and getters
		void set_Z0(const double& Z0) { m_Z0 = Z0; };
		void set_er(const double& er) { m_er = er; };
		void set_d(const double& d) { m_d = d; };
		void set_Cn(const std::vector<double>& Cn) { m_Cn = Cn; };

		double get_Z0() const { return m_Z0; };
		double get_er() const { return m_er; };
		double get_d() const { return m_d; };
		std::vector<double> get_Cn() const { return m_Cn; };
	};

}