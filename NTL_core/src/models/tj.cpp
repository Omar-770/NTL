#include "TJ.h"

namespace TJ
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
        if (std::abs(B) < 1e-12)
        {
            // If exactly 0, default to positive real, else maintain the phase
            B = (B == std::complex<double>(0, 0)) ? 1e-12 : (B / std::abs(B)) * 1e-12;
        }

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


    matrix3x3cd calculate_Y_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2, int M2,
        double d3, const std::vector<double>& Cn3, int M3, double f, int K)
    {
        matrix2x2cd y2, y3;
        y2 = NTL::calculate_Y_matrix(Z0, er, d2, Cn2, M2, f, K);
        y3 = NTL::calculate_Y_matrix(Z0, er, d3, Cn3, M3, f, K);

        matrix3x3cd Y = matrix3x3cd::Zero();

        Y(0, 0) += y2(0, 0);
        Y(0, 1) += y2(0, 1);
        Y(1, 0) += y2(1, 0);
        Y(1, 1) += y2(1, 1);

        Y(0, 0) += y3(0, 0);
        Y(0, 2) += y3(0, 1);
        Y(2, 0) += y3(1, 0);
        Y(2, 2) += y3(1, 1);

        return Y;
    }

    matrix3x3cd calculate_Y_matrix(const TJ& tj, double f, int K)
    {
        return calculate_Y_matrix(tj.get_Z0(), tj.get_er(), tj.get_ntl2().get_d(),
            tj.get_ntl2().get_Cn(), tj.get_ntl2().get_M(), tj.get_ntl3().get_d(), tj.get_ntl3().get_Cn(),
            tj.get_ntl2().get_M(), f, K);
    }


    matrix3x3cd calculate_S_matrix(double Z0, double er, double d2, const std::vector<double>& Cn2, int M2,
        double d3, const std::vector<double>& Cn3, int M3, double f, std::array<std::complex<double>, 3> Zl, int K)
    {
        matrix3x3cd Y_wpd = calculate_Y_matrix(Z0, er, d2, Cn2, M2, d3, Cn3, M3, f, K);

        matrix3x3cd Y_ref = matrix3x3cd::Zero();
        matrix3x3cd Y_ref_conj = matrix3x3cd::Zero();
        matrix3x3cd F = matrix3x3cd::Zero();
        matrix3x3cd F_inv = matrix3x3cd::Zero();

        for (int i = 0; i < 3; ++i)
        {
            std::complex<double> y_i = 1.0 / Zl[i];
            Y_ref(i, i) = y_i;
            Y_ref_conj(i, i) = std::conj(y_i);

            // Power waves require the real part of the reference admittance
            double g_i = std::real(y_i);

            // Safety clamp to prevent division by zero or sqrt of negative numbers
            // if the termination is purely reactive (which shouldn't happen physically, but protects the math)
            if (g_i < 1e-15) g_i = 1e-15;

            F(i, i) = 1.0 / std::sqrt(g_i);
            F_inv(i, i) = std::sqrt(g_i);
        }

        // Full Kurokawa N-port transformation
        matrix3x3cd S = F * (Y_ref_conj - Y_wpd) * (Y_ref + Y_wpd).inverse() * F_inv;

        return S;
    }


    matrix3x3cd calculate_S_matrix(const TJ& tj, double f, std::array<std::complex<double>, 3> Zl, int K)
    {
        return calculate_S_matrix(tj.get_Z0(), tj.get_er(), tj.get_ntl2().get_d(),
            tj.get_ntl2().get_Cn(), tj.get_ntl2().get_M(), tj.get_ntl3().get_d(), tj.get_ntl3().get_Cn(),
            tj.get_ntl3().get_M(), f, Zl, K);

    }


    std::pair<matrix3x3cd, std::vector<matrix3x3cd>> calculate_S_matrix_with_grad(
        double Z0, double er, double d2, const std::vector<double>& Cn2, int M2,
        double d3, const std::vector<double>& Cn3, int M3, double f,
        std::array<std::complex<double>, 3> Zl, int K)
    {
        // 1. Calculate NTL derivatives
        auto [T2, dT2] = NTL::calculate_T_matrix_with_grad(Z0, er, d2, Cn2, M2, f, K);
        auto [T3, dT3] = NTL::calculate_T_matrix_with_grad(Z0, er, d3, Cn3, M3, f, K);

        // 2. Convert to Y derivatives
        auto [y2, dy2] = T_to_Y_with_grad(T2, dT2);
        auto [y3, dy3] = T_to_Y_with_grad(T3, dT3);

        // 3. Construct Global Y Matrix

        matrix3x3cd Y = matrix3x3cd::Zero();

        // Port 1 (Node 0) connects to Port 1 of both arms (y2(0,0), y3(0,0))
        Y(0, 0) = y2(0, 0);
        Y(0, 1) = y2(0, 1);
        Y(0, 2) = y3(0, 1);

        // Port 2 (Node 1) connects to Arm 2 out and Resistor
        Y(1, 0) = y2(1, 0);
        Y(1, 1) = y2(1, 1);

        // Port 3 (Node 2) connects to Arm 3 out and Resistor
        Y(2, 0) = y3(1, 0);
        Y(2, 2) = y3(1, 1);

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
        // 5. Calculate S Matrix and Gradient (Full Kurokawa Formulation)
        matrix3x3cd Y_ref = matrix3x3cd::Zero();
        matrix3x3cd Y_ref_conj = matrix3x3cd::Zero();
        matrix3x3cd F = matrix3x3cd::Zero();
        matrix3x3cd F_inv = matrix3x3cd::Zero();

        for (int i = 0; i < 3; ++i)
        {
            std::complex<double> y_i = 1.0 / Zl[i];
            Y_ref(i, i) = y_i;
            Y_ref_conj(i, i) = std::conj(y_i);

            double g_i = std::real(y_i);
            if (g_i < 1e-15) g_i = 1e-15;

            F(i, i) = 1.0 / std::sqrt(g_i);
            F_inv(i, i) = std::sqrt(g_i);
        }

        matrix3x3cd _M = Y_ref + Y; // Term (Yref + Y)
        matrix3x3cd M_inv = _M.inverse();

        // Exact Kurokawa Power S-Matrix
        matrix3x3cd S = F * (Y_ref_conj - Y) * M_inv * F_inv;

        matrix3x3cd I = matrix3x3cd::Identity();

        // The analytical Kurokawa Jacobian Prefactor
        matrix3x3cd PreFactor = -(I + S) * F;

        std::vector<matrix3x3cd> dS(total_params);

        for (size_t i = 0; i < total_params; ++i)
        {
            dS[i] = PreFactor * dY[i] * M_inv * F_inv;
        }

        return { S, dS };
    }

    /// CLASS_TJ_BEGIN

    matrix3x3cd TJ::Y_matrix(double f, int K) const
    {
        return calculate_Y_matrix(*this, f, K);
    }
    matrix3x3cd TJ::S_matrix(double f, std::array<std::complex<double>, 3> Zl, int K) const
    {
        return calculate_S_matrix(*this, f, Zl, K);
    }

}