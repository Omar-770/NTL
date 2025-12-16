#include "optimisation/WPD_opt_2.h"

namespace WPD
{
    opt_2::opt_2(const opt_setup& setup) : optimiser(setup)
    {
        m_Z0 = setup.Z0;
        m_er = setup.er;
        m_d = setup.d;
        m_M = setup.M;
        m_Zref = setup.Zref;
        m_freqs = setup.freqs;
        m_K = setup.K;
        m_Z_min = setup.Z_min;
        m_Z_max = setup.Z_max;
        m_R_min = setup.R_min;
        m_R_max = setup.R_max;
        m_split = setup.split;

        if (m_M > m_N)
            throw(std::invalid_argument("Sine terms must be fewer or equal to the number of terms"));
        if (m_split.size() < m_freqs.size() && m_split.size() != 1)
            throw(std::invalid_argument("Number of splits & frequency points mismatch"));
        if (m_split.size() == 1)
        {
            m_split.resize(m_freqs.size());
            for (int i = 1; i < m_split.size(); i++)
                m_split[i] = m_split[0];
        }
        if (m_Z0 < 1e-6 || m_er < 1e-6 || m_d < 1e-6)
            throw(std::invalid_argument("Invalid NTL physical charachteristics"));
        if (m_freqs.empty())
            throw(std::invalid_argument("Empty frequency vector"));
        if (m_Z_min < 1e-6 || m_Z_max < 1e-6 || m_Z_max < m_Z_min)
            throw(std::invalid_argument("Invalid min/max impedance(s)"));
        if (m_R_min < 1e-6 || m_R_max < 1e-6 || m_R_max < m_R_min)
            throw(std::invalid_argument("Invalid min/max resistance"));

        m_F = m_freqs.size();
        m_N1 = m_N;

        auto set_physical_params = [&](NTL::NTL& ntl)
            {
                ntl.set_Z0(m_Z0);
                ntl.set_er(m_er);
                ntl.set_d(m_d);
                ntl.set_M(m_M);
            };

        for (auto& ntl : m_ntl)
            set_physical_params(ntl);

        m_wpd.set_Z0(m_Z0);
        m_wpd.set_er(m_er);
        m_wpd.set_M(m_M);

        auto resize_vec = [&](std::vector<std::complex<double>>& vec)
            {
                vec.resize(m_F);
            };

        resize_vec(m_arm2_Zl);
        resize_vec(m_arm2_Zs);
        resize_vec(m_arm3_Zl);
        resize_vec(m_arm3_Zs);

        m_d_arms = false;
        m_d_outputs = false;

        int hw_threads = std::thread::hardware_concurrency();
        if (hw_threads == 0) hw_threads = 4; 

        omp_set_num_threads(std::min<int>(m_F, hw_threads));
    }

    opt_result opt_2::optimise(console mode)
    {
        m_out = mode == console::active ? true : false;
        m_use_local_grad_free = true;

        if (m_out)
        {
            std::cout << "Initialising WPD optimisation...\n";
            std::cout << "Optimising arms jointly...\n\n";
        }

        optimise_arms();

        if (m_out)
            std::cout << "\n\nOptimising output transformers...\n\n";

        optimise_transformers();

        optimise_R();

        m_wpd.set_arms(m_ntl[0], m_ntl[1]);
        m_wpd.set_R(m_R);

        return opt_result{ {}, m_wpd, m_ntl[2], m_ntl[3] };
    }

    void opt_2::print_result_logs()
    {

        for (int i = 0; i < m_F; i++)
        {
            double f = m_freqs[i];
            std::cout << "[freq: " << f << " ]\n";
            std::cout << "NTL2: (theta) = " << m_ntl[0].electrical_length(f) << '\n';
            std::cout << "NTL3: (theta) = " << m_ntl[1].electrical_length(f) << '\n';
            std::cout << "NTL2: (H) = " << m_ntl[0].V_transfer(f, m_arm2_Zl[i]) << '\n';
            std::cout << "NTL3: (H) = " << m_ntl[1].V_transfer(f, m_arm3_Zl[i]) << '\n';
            std::cout << "NTL2: (Gamma) = " << m_ntl[0].S11(f, m_arm2_Zs[i], m_arm2_Zl[i]) << '\n';
            std::cout << "NTL3: (Gamma) = " << m_ntl[1].S11(f, m_arm3_Zs[i], m_arm3_Zl[i]) << '\n';
            std::cout << "Leaked power (Optimised R) = " <<
                10 * std::log10(std::pow(std::abs(m_ntl[0].V_transfer(f, m_Zref) - m_ntl[1].V_transfer(f, m_Zref)), 2) / (2 * m_R))
                << " dBm\n\n";

        }

        std::cout << "m_arm2_Zs: " << m_arm2_Zs << std::endl;
        std::cout << "m_arm2_Zl: " << m_arm2_Zl << std::endl;
        std::cout << "m_arm3_Zs: " << m_arm3_Zs << std::endl;
        std::cout << "m_arm3_Zl: " << m_arm3_Zl << std::endl;

        std::cout << "Final WPD:\n";
        std::cout << "Arm2: (Cn) = " << m_ntl[0].get_Cn() << '\n';
        std::cout << "Arm3: (Cn) = " << m_ntl[1].get_Cn() << '\n';
        std::cout << "Out2: (Cn) = " << m_ntl[2].get_Cn() << '\n';
        std::cout << "Out3: (Cn) = " << m_ntl[3].get_Cn() << '\n';
        std::cout << "R: = " << m_R << '\n';
    }

    void opt_2::optimise_arms()
    {
        m_N *= 2;
        m_lb.assign(m_N, -1);
        m_ub.assign(m_N, 1);
        m_toll_z.assign(m_K, 1e-6);
        m_toll_bounds.assign(4, 1e-6);

        for (int i = 0; i < m_F; i++)
        {
            m_arm2_Zs[i] = m_Zref * (1 + m_split[i]);
            m_arm2_Zl[i] = m_Zref * std::sqrt(m_split[i]);
            m_arm3_Zs[i] = m_Zref * (1 + 1.0 / m_split[i]);
            m_arm3_Zl[i] = m_Zref / std::sqrt(m_split[i]);
        }

        std::vector<double> Cn;

        if (m_d_arms)
        {
            double init_d = m_d;
            int init_max_attempts = m_max_attempts;
            double resolution = 1e-3;

            if (m_out)
                std::cout << "\n\t===>>> Starting with d = " << m_d * 1000 << "mm" << std::endl;
            optimiser_result init_result = run_optimiser(m_out ? console::active : console::inactive);
            optimiser_result new_result = init_result;

            if (init_result.final_error < m_accepted_error)
            {
                double error_this_attempt = init_result.final_error;

                while (true)
                {
                    m_d -= resolution;

                    if (m_out)
                        std::cout << "\n\t===>>> Trimming NTL by " << resolution * 1000 << "mm to:\t"
                        << m_d * 1000 << "mm" << std::endl;

                    optimiser_result result_this_attempt = run_optimiser(m_out ? console::active : console::inactive);

                    if (result_this_attempt.final_error > m_accepted_error)
                    {
                        m_d += resolution;
                        break;
                    }

                    new_result = result_this_attempt;
                }

            }
            m_ntl[0].set_d(m_d);
            m_ntl[1].set_d(m_d);
            m_d = init_d;
            m_max_attempts = init_max_attempts;
        }
        else
        {
            Cn = run_optimiser(m_out ? console::active : console::inactive).optimised_cn;
        }

        m_ntl[0].set_Cn(std::vector<double>(Cn.begin(), Cn.begin() + m_N1));
        m_ntl[1].set_Cn(std::vector<double>(Cn.begin() + m_N1, Cn.begin() + 2 * m_N1));

        if (m_out && m_d_arms)
            std::cout << "Final NTL length:\t" << m_ntl[0].get_d() * 1000 << "mm" << std::endl;


    }

    void opt_2::optimise_transformers()
    {
        //creates two NTL::opt objects to make the transformers
        NTL::opt_setup setup;
        setup.N = m_N1;
        setup.lb.assign(m_N1, -1);
        setup.ub.assign(m_N1, 1);
        setup.toll_bounds.assign(2, 1e-6);
        setup.toll_z.assign(m_K, 1e-6);
        setup.GBL_MAX = m_GBL_MAX;
        setup.LCL_MAX = m_LCL_MAX;
        setup.accepted_error = m_accepted_error;
        setup.max_attempts = m_max_attempts;
        setup.Z0 = m_Z0;
        setup.er = m_er;
        setup.d = m_d;
        setup.M = m_M;
        setup.freqs = m_freqs;
        setup.K = m_K;
        setup.Z_min = m_Z_min;
        setup.Z_max = m_Z_max;
        setup.Z_at_0 = m_Zref;
        setup.Z_at_d = m_Zref;

        if (m_out) std::cout << "\n\t-->Optimising first transformer...\n\n";

        setup.Zs.resize(m_F);
        setup.Zl.resize(m_F);

        for (int i = 0; i < m_F; i++) {
            setup.Zs[i] = m_Zref * std::sqrt(m_split[i]);
            setup.Zl[i] = m_Zref;
        }


        NTL::opt opt2(setup);
        if (m_d_outputs)
            m_ntl[2] = opt2.optimise_d(m_out ? console::active : console::inactive).ntl;
        else
            m_ntl[2] = opt2.optimise(m_out ? console::active : console::inactive).ntl;

        if (m_out) std::cout << "\n\t-->Optimising second transformer...\n\n";
        for (int i = 0; i < m_F; i++) {
            setup.Zs[i] = m_Zref / std::sqrt(m_split[i]);
            setup.Zl[i] = m_Zref;
        }


        NTL::opt opt3(setup);
        if (m_d_outputs)
            m_ntl[3] = opt3.optimise_d(m_out ? console::active : console::inactive).ntl;
        else
            m_ntl[3] = opt3.optimise(m_out ? console::active : console::inactive).ntl;
    }

    void opt_2::optimise_R()
    {
        // Save state
        bool original_out = m_out;
        m_out = false;

        if (original_out)
            std::cout << "Calculating optimal isolation resistance (Hybrid Sweep + Golden Section)...\n";

        // ---------------------------------------------------------
        // 1. Pre-calculate the Series Impedance for all frequencies
        // ---------------------------------------------------------
        std::vector<std::complex<double>> Z_series_vec;
        Z_series_vec.reserve(m_F);

        for (int i = 0; i < m_F; i++)
        {
            double f = m_freqs[i];

            // Impedance looking into the Transformers (Load = Zref)
            std::complex<double> Z_trans2 = m_ntl[2].Zin(m_Zref, f, m_K);
            std::complex<double> Z_trans3 = m_ntl[3].Zin(m_Zref, f, m_K);

            // Impedance looking back into the Arms (Source = Short/Ground for Odd Mode)
            std::complex<double> Z_arm2_back = m_ntl[0].Zout(0.0, f, m_K);
            std::complex<double> Z_arm3_back = m_ntl[1].Zout(0.0, f, m_K);

            // Parallel combination at each node (Transformer || Arm)
            std::complex<double> Y_node1 = (1.0 / Z_trans2) + (1.0 / Z_arm2_back);
            std::complex<double> Y_node2 = (1.0 / Z_trans3) + (1.0 / Z_arm3_back);

            // Total Series Impedance seen by the resistor (Node1 + Node2)
            std::complex<double> Z_series = (1.0 / Y_node1) + (1.0 / Y_node2);
            Z_series_vec.push_back(Z_series);
        }

        // Cost function: Sum of squared reflection coefficients across the band
        auto calculate_cost = [&](double R) -> double {
            double cost = 0.0;
            for (const auto& Z : Z_series_vec)
            {
                // |Gamma|^2 = |(Z - R)/(Z + R)|^2
                cost += std::norm((Z - R) / (Z + R));
            }
            return cost;
            };

        // ---------------------------------------------------------
        // 2. Coarse Sweep (Find Global Minimum Region)
        // ---------------------------------------------------------
        double best_R_guess = m_R_min;
        double min_cost = std::numeric_limits<double>::max();
        int sweep_points = 100;
        double step = (m_R_max - m_R_min) / sweep_points;

        for (int i = 0; i <= sweep_points; ++i)
        {
            double r = m_R_min + i * step;
            double cost = calculate_cost(r);
            if (cost < min_cost)
            {
                min_cost = cost;
                best_R_guess = r;
            }
        }

        // ---------------------------------------------------------
        // 3. Golden Section Search (Refine for Precision)
        // ---------------------------------------------------------
        // Bracket the search around the best guess
        double bracket_width = step * 2.0;
        double a = std::max(m_R_min, best_R_guess - bracket_width);
        double b = std::min(m_R_max, best_R_guess + bracket_width);

        double gr = (std::sqrt(5.0) + 1.0) / 2.0;
        double c = b - (b - a) / gr;
        double d = a + (b - a) / gr;

        while (std::abs(b - a) > 1e-4) // Precision tolerance
        {
            if (calculate_cost(c) < calculate_cost(d)) {
                b = d;
                d = c;
                c = b - (b - a) / gr;
            }
            else {
                a = c;
                c = d;
                d = a + (b - a) / gr;
            }
        }

        m_R = (a + b) / 2.0;

        // Restore state
        m_out = original_out;
        if (m_out)
            std::cout << "  >> Final Optimal R: " << m_R << " Ohms\n";
    }

    double opt_2::min_objective(const std::vector<double>& Cn) const
    {
        double sum_squares{};
        std::vector<double> Cn2(Cn.begin(), Cn.begin() + m_N1);
        std::vector<double> Cn3(Cn.begin() + m_N1, Cn.begin() + 2 * m_N1);

#pragma omp parallel        
        {
            double thread_sum_squares{};
#pragma omp for nowait
            for (int i = 0; i < m_F; i++)
            {
                double f = m_freqs[i];
                std::complex<double> Zl2 = m_arm2_Zl[i];
                std::complex<double> Zs2 = m_arm2_Zs[i];
                std::complex<double> Zl3 = m_arm3_Zl[i];
                std::complex<double> Zs3 = m_arm3_Zs[i];

                matrix2x2cd T2 = NTL::calculate_T_matrix(m_Z0, m_er, m_d, Cn2, m_M, f, m_K);
                std::complex<double> A2 = T2(0, 0), B2 = T2(0, 1), C2 = T2(1, 0), D2 = T2(1, 1);
                std::complex<double> Zin2 = (A2 * Zl2 + B2) / (C2 * Zl2 + D2);
                std::complex<double> H2 = 1.0 / (A2 + B2 / Zl2);
                std::complex<double> gamma2 = (Zs2 - std::conj(Zin2)) / (Zs2 + Zin2);

                matrix2x2cd T3 = NTL::calculate_T_matrix(m_Z0, m_er, m_d, Cn3, m_M, f, m_K);
                std::complex<double> A3 = T3(0, 0), B3 = T3(0, 1), C3 = T3(1, 0), D3 = T3(1, 1);
                std::complex<double> Zin3 = (A3 * Zl3 + B3) / (C3 * Zl3 + D3);
                std::complex<double> H3 = 1.0 / (A3 + B3 / Zl3);
                std::complex<double> gamma3 = (Zs3 - std::conj(Zin3)) / (Zs3 + Zin3);

                thread_sum_squares += std::pow(std::norm(gamma2), 2);
                thread_sum_squares += std::pow(std::norm(gamma3), 2);
                thread_sum_squares += std::pow(H2.real() - H3.real(), 2);
                thread_sum_squares += std::pow(H2.imag() - H3.imag(), 2);
            }

#pragma omp critical
            {
                sum_squares += thread_sum_squares;
            }
        }

        return sum_squares;
    }

    void opt_2::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
    {
        //terminal impedances = Zref
        const double* Cn2 = Cn;
        const double* Cn3 = Cn + m_N1;

        res[0] = NTL::calculate_Z(Cn2, m_N1, m_M, m_Z0, m_d, 0.0) - m_Zref;
        res[1] = NTL::calculate_Z(Cn2, m_N1, m_M, m_Z0, m_d, m_d) - m_Zref;
        res[2] = NTL::calculate_Z(Cn3, m_N1, m_M, m_Z0, m_d, 0.0) - m_Zref;
        res[3] = NTL::calculate_Z(Cn3, m_N1, m_M, m_Z0, m_d, m_d) - m_Zref;
    }

    void opt_2::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
    {
        //All K subsections impedances < Z_max
        const double* Cn2 = Cn;
        const double* Cn3 = Cn + m_N1;

        int half_K = m_K / 2;
        double dz = m_d / m_K;

        for (int i = 0; i < half_K; i++)
            res[i] = NTL::calculate_Z(Cn2, m_N1, m_M, m_Z0, m_d, (double(i) + 0.5) * dz) - m_Z_max;

        for (int i = 0; i < half_K; i++)
            res[half_K + i] = NTL::calculate_Z(Cn3, m_N1, m_M, m_Z0, m_d, (double(i) + 0.5) * dz) - m_Z_max;
    }

    void opt_2::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
    {
        //All K subsections impedances > Z_min
        const double* Cn2 = Cn;
        const double* Cn3 = Cn + m_N1;

        int half_K = m_K / 2;
        double dz = m_d / m_K;

        for (int i = 0; i < half_K; i++)
            res[i] = m_Z_min - NTL::calculate_Z(Cn2, m_N1, m_M, m_Z0, m_d, (double(i) + 0.5) * dz);

        for (int i = 0; i < half_K; i++)
            res[half_K + i] = m_Z_min - NTL::calculate_Z(Cn3, m_N1, m_M, m_Z0, m_d, (double(i) + 0.5) * dz);
    }

    double opt_2::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
    {
        //No need m_use_local_grad_free is true
        return 0.0;
    }
}