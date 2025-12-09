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
	struct DATA;

	using matrix2x2cd = Eigen::Matrix<std::complex<double>, 2, 2>;

	inline constexpr double M_PI = 3.14159265358979311599796346854;
	inline constexpr double M_C = 299792458;

	double calculate_Z(double Z0, double d, const std::vector<double>& Cn, int M, double z);
	double calculate_Z(const NTL& ntl, double z);
	double calculate_Z(const double* Cn, const size_t& n, const size_t& m, const double& Z0, const double& d, const double& z); //for nlopt

	std::complex<double> calculate_Zin(double Z0, double e_r, double d, const std::vector<double>& Cn,
		int M, double Zl, double f, int K = 50);
	std::complex<double> calculate_Zin(const NTL& ntl, double Zl, double f, int K = 50);

	double calculate_W_H(double z, double e_r);
	double calculate_er_eff(double z, double e_r);

	matrix2x2cd calculate_T_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, int M,
		double f, int K = 50);
	matrix2x2cd calculate_T_matrix(const NTL& ntl, double f, int K = 50);

	std::pair<matrix2x2cd, std::vector<matrix2x2cd>> calculate_T_matrix_with_grad(
		double Z0, double er, double d, const std::vector<double>& Cn, int M, double f, int K = 50);
	

	matrix2x2cd calculate_S_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, int M,
		double f, double Zs, double Zl, int K = 50);
	matrix2x2cd calculate_S_matrix(const NTL& ntl, double f, double Zs, double Zl, int K = 50);

	matrix2x2cd calculate_Y_matrix(double Z0, double e_r, double d, const std::vector<double>& Cn, int M,
		double f, int K = 50);
	matrix2x2cd calculate_Y_matrix(const NTL& ntl, double f, int K = 50);


	struct DATA
	{
		double Z0;
		double er;
		double d;

		std::vector<double> Cn;
		int M;
	};

	class NTL
	{
	public:
		NTL() : m_Z0(0), m_er(0), m_d(0), m_Cn({ 0 }), m_M(0) {

		}

		NTL(double Z0, double er, double d)
			: m_Z0(Z0), m_er(er), m_d(d) {
			m_Cn = std::vector<double>{ 0.0 };
			m_M = 0;
		};

		NTL(double Z0, double er, double d, const std::vector<double>& Cn, int M = 0)
			: m_Z0(Z0), m_er(er), m_d(d), m_Cn(Cn), m_M(M) {
		};


		NTL(const DATA& data)
			: m_Z0(data.Z0), m_er(data.er), m_d(data.d), m_Cn(data.Cn), m_M(data.M) {
		};

		//impedance profiles
		double Z(double z) const;
		std::complex<double> Zin(double Zl, double f, int K = 50) const;

		double W_H(double z) const;
		double er_eff(const double& z) const;

		//data extraction
		std::vector<std::pair<double, double>> get_Z_vec(double step_size = 0.1) const;
		std::vector<std::pair<double, double>> get_w_h_vec(double step_size = 0.1) const;

		DATA data() const;

		//microwave matrices
		matrix2x2cd T_matrix(double f, int K = 50) const;
		matrix2x2cd S_matrix(double f, double Zs, double Zl, int K = 50) const;
		matrix2x2cd Y_matrix(double f, int K = 50) const;

		//reflection coefficient
		std::complex<double> S11(double f, double Zs, double Zl, int K = 50) const;

	private:

		//Physical parameters
		double m_Z0; // Characteristic impedance
		double m_er; // Relative permittivity
		double m_d;  // Length of the line

		std::vector<double> m_Cn; // Coefficients vector
		int m_M; // Number of sine terms

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
		int get_M() const { return m_M; };
	};

}