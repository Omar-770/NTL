#include "wpd.h"

namespace WPD
{
    // Helper: Convert T-matrix and its gradients to Y-matrix and its gradients
    // Y = [ D/B,  -1/B ]
    //     [ -1/B,  A/B ] 
    // assuming AD-BC = 1 (reciprocal/lossless)
    static std::pair<matrix2x2cd, std::vector<matrix2x2cd>> T_to_Y_with_grad(
        const matrix2x2cd& T, const std::vector<matrix2x2cd>& dT)
    {
        std::complex<double> A = T(0, 0);
        std::complex<double> B = T(0, 1);
        std::complex<double> D = T(1, 1);

        // Guard against division by zero (B=0)
        if (std::abs(B) < 1e-12) B = 1e-12; // Simple numerical stability patch

        std::complex<double> invB = 1.0 / B;
        std::complex<double> invB2 = invB * invB;

        matrix2x2cd Y;
        Y(0, 0) = D * invB;
        Y(0, 1) = -invB;
        Y(1, 0) = -invB;
        Y(1, 1) = A * invB;

        size_t N = dT.size();
        std::vector<matrix2x2cd> dY(N);

        for (size_t i = 0; i < N; ++i)
        {
            std::complex<double> dA = dT[i](0, 0);
            std::complex<double> dB = dT[i](0, 1);
            std::complex<double> dD = dT[i](1, 1);

            // Quotient rule derivatives
            // d(1/B) = -dB / B^2
            std::complex<double> d_invB = -dB * invB2;

            // d(D/B) = (dD*B - D*dB) / B^2
            std::complex<double> d_DB = (dD * B - D * dB) * invB2;

            // d(A/B) = (dA*B - A*dB) / B^2
            std::complex<double> d_AB = (dA * B - A * dB) * invB2;

            dY[i](0, 0) = d_DB;
            dY[i](0, 1) = d_invB;
            dY[i](1, 0) = d_invB;
            dY[i](1, 1) = d_AB;
        }

        return { Y, dY };
    }

    std::pair<matrix3x3cd, std::vector<matrix3x3cd>> calculate_S_matrix_with_grad(
        double Z0, double er, double d2, const std::vector<double>& Cn2,
        double d3, const std::vector<double>& Cn3, int M, double R, double f,
        std::array<std::complex<double>, 3> Zl, int K)
    {
        // 1. Calculate NTL derivatives
        auto [T2, dT2] = NTL::calculate_T_matrix_with_grad(Z0, er, d2, Cn2, M, f, K);
        auto [T3, dT3] = NTL::calculate_T_matrix_with_grad(Z0, er, d3, Cn3, M, f, K);

        // 2. Convert to Y derivatives
        auto [y2, dy2] = T_to_Y_with_grad(T2, dT2);
        auto [y3, dy3] = T_to_Y_with_grad(T3, dT3);

        // 3. Construct Global Y Matrix
        // Resistor conductance G = 1/R
        double G = 1.0 / R;
        double dG_dR = -1.0 / (R * R);

        matrix2x2cd y_R;
        y_R << G, -G, -G, G;

        matrix3x3cd Y = matrix3x3cd::Zero();

        // Port 1 (Node 0) connects to Port 1 of both arms (y2(0,0), y3(0,0))
        Y(0, 0) = y2(0, 0) + y3(0, 0);
        Y(0, 1) = y2(0, 1);
        Y(0, 2) = y3(0, 1);

        // Port 2 (Node 1) connects to Arm 2 out and Resistor
        Y(1, 0) = y2(1, 0);
        Y(1, 1) = y2(1, 1) + y_R(0, 0);
        Y(1, 2) = y_R(0, 1);

        // Port 3 (Node 2) connects to Arm 3 out and Resistor
        Y(2, 0) = y3(1, 0);
        Y(2, 1) = y_R(1, 0);
        Y(2, 2) = y3(1, 1) + y_R(1, 1);

        // 4. Construct Global Gradient Vector
        // Structure: [ ...Cn2... | ...Cn3... | R ]
        size_t n2_size = Cn2.size();
        size_t n3_size = Cn3.size();
        size_t total_params = n2_size + n3_size + 1;

        std::vector<matrix3x3cd> dY(total_params, matrix3x3cd::Zero());

        // Fill dY for Cn2
        for (size_t i = 0; i < n2_size; ++i)
        {
            dY[i](0, 0) += dy2[i](0, 0);
            dY[i](0, 1) += dy2[i](0, 1);
            dY[i](1, 0) += dy2[i](1, 0);
            dY[i](1, 1) += dy2[i](1, 1);
        }

        // Fill dY for Cn3
        for (size_t i = 0; i < n3_size; ++i)
        {
            size_t idx = n2_size + i;
            dY[idx](0, 0) += dy3[i](0, 0);
            dY[idx](0, 2) += dy3[i](0, 1);
            dY[idx](2, 0) += dy3[i](1, 0);
            dY[idx](2, 2) += dy3[i](1, 1);
        }

        // Fill dY for R (last element)
        size_t r_idx = total_params - 1;
        dY[r_idx](1, 1) += dG_dR;
        dY[r_idx](1, 2) -= dG_dR;
        dY[r_idx](2, 1) -= dG_dR;
        dY[r_idx](2, 2) += dG_dR;

        // 5. Calculate S Matrix and Gradient
        // S = (Yref - Y)(Yref + Y)^-1
        // dS = -(I + S) * dY * (Yref + Y)^-1

        matrix3x3cd Z_ref = matrix3x3cd::Zero();
        Z_ref(0, 0) = Zl[0];
        Z_ref(1, 1) = Zl[1];
        Z_ref(2, 2) = Zl[2];
        matrix3x3cd Y_ref = Z_ref.inverse();

        matrix3x3cd _M = Y_ref + Y; // Term (Yref + Y)
        matrix3x3cd M_inv = _M.inverse();

        matrix3x3cd S = (Y_ref - Y) * M_inv;
        matrix3x3cd I = matrix3x3cd::Identity();
        matrix3x3cd PreFactor = -(I + S);

        std::vector<matrix3x3cd> dS(total_params);

        for (size_t i = 0; i < total_params; ++i)
        {
            dS[i] = PreFactor * dY[i] * M_inv;
        }

        return { S, dS };
    }
}

namespace WPD
{
	matrix3x3cd calculate_Y_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2,
		double d3, const std::vector<double>& Cn3, int M, double R, double f, int K)
	{
		double G = 1 / R;

		matrix2x2cd y2, y3;
		y2 = NTL::calculate_Y_matrix(Z0, er, d2, Cn2, M, f, K);
		y3 = NTL::calculate_Y_matrix(Z0, er, d3, Cn3, M, f, K);
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
			wpd.get_ntl2().get_Cn(), wpd.get_ntl3().get_d(), wpd.get_ntl3().get_Cn(), wpd.get_M(),
			wpd.get_R(), f, K);
	}

	matrix3x3cd calculate_S_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2,
		double d3, const std::vector<double>& Cn3, int M, double R, double f, std::array<std::complex<double>, 3> Zl, int K)
	{
		matrix3x3cd Y_wpd = calculate_Y_matrix(Z0, er, d2, Cn2, d3, Cn3, M, R, f, K);
		matrix3x3cd Z_ref = matrix3x3cd::Zero();
		Z_ref(0, 0) = Zl[0];
		Z_ref(1, 1) = Zl[1];
		Z_ref(2, 2) = Zl[2];

		matrix3x3cd Y_ref = Z_ref.inverse();
		
		matrix3x3cd S = (Y_ref - Y_wpd) * (Y_ref + Y_wpd).inverse();


		return S;
	}

	matrix3x3cd calculate_S_matrix(const WPD& wpd, double f, std::array<std::complex<double>, 3> Zl, int K)
	{
		return calculate_S_matrix(wpd.get_Z0(), wpd.get_er(), wpd.get_ntl2().get_d(),
			wpd.get_ntl2().get_Cn(), wpd.get_ntl3().get_d(), wpd.get_ntl3().get_Cn(), wpd.get_M(),
			wpd.get_R(), f, Zl, K);
	}

	/// CLASS_WPD_BEGIN

	matrix3x3cd WPD::Y_matrix(double f, int K) const
	{
		return calculate_Y_matrix(*this, f, K);
	}
	matrix3x3cd WPD::S_matrix(double f, std::array<std::complex<double>, 3> Zl, int K) const
	{
		return calculate_S_matrix(*this, f, Zl, K);
	}
}