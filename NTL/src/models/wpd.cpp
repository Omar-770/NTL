#include "wpd.h"

namespace NTL
{
	matrix3x3cd calculate_Y_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2,
		double d3, const std::vector<double>& Cn3, double R, double f, int K)
	{
		double G = 1 / R;

		matrix2x2cd y2, y3;
		y2 = calculate_Y_matrix(Z0, er, d2, Cn2, f, K);


		y3 = calculate_Y_matrix(Z0, er, d3, Cn3, f, K);
		matrix2x2cd y_R
		{
			{G, -G},
			{-G, G}
		};


		matrix3x3cd Y = matrix3x3cd::Zero();

		Y(0, 0) += y2(0, 0);
		Y(0, 1) += y2(0, 1);
		Y(1, 0) += y2(1, 0);
		Y(1, 1) += y2(1, 1);

		Y(0, 0) += y3(0, 0);
		Y(0, 2) += y3(0, 1);
		Y(2, 0) += y3(1, 0);
		Y(2, 2) += y3(1, 1);

		Y(1, 1) += y_R(0, 0);
		Y(1, 2) += y_R(0, 1);
		Y(2, 1) += y_R(1, 0);
		Y(2, 2) += y_R(1, 1);

		return Y;
	}

	matrix3x3cd calculate_Y_matrix(const WPD& wpd, double f, int K)
	{
		return calculate_Y_matrix(wpd.get_Z0(), wpd.get_er(), wpd.get_ntl2().get_d(),
			wpd.get_ntl2().get_Cn(), wpd.get_ntl3().get_d(), wpd.get_ntl3().get_Cn(),
			wpd.get_R(), f, K);
	}

	matrix3x3cd calculate_S_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2,
		double d3, const std::vector<double>& Cn3, double R, double f, double Z1, double Z2, double Z3, int K)
	{
		matrix3x3cd Y_wpd = calculate_Y_matrix(Z0, er, d2, Cn2, d3, Cn3, R, f, K);

		matrix3x3cd Y_ref = matrix3x3cd::Zero();
		Y_ref(0, 0) = 1 / Z1; Y_ref(1, 1) = 1 / Z2; Y_ref(2, 2) = 1 / Z3;

		matrix3x3cd S = (Y_ref - Y_wpd) * (Y_ref + Y_wpd).inverse();


		return S;
	}

	matrix3x3cd calculate_S_matrix(const WPD& wpd, double f, double Z1, double Z2, double Z3, int K)
	{
		return calculate_S_matrix(wpd.get_Z0(), wpd.get_er(), wpd.get_ntl2().get_d(),
			wpd.get_ntl2().get_Cn(), wpd.get_ntl3().get_d(), wpd.get_ntl3().get_Cn(),
			wpd.get_R(), f, Z1, Z2, Z3, K);
	}

	/// CLASS_WPD_BEGIN

	matrix3x3cd WPD::Y_matrix(double f, int K) const
	{
		return calculate_Y_matrix(*this, f, K);
	}
	matrix3x3cd WPD::S_matrix(double f, double Z1, double Z2, double Z3, int K) const
	{
		return calculate_S_matrix(*this, f, Z1, Z2, Z3, K);
	}
}