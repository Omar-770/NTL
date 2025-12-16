#include "WPD_opt.h"

namespace WPD
{
    opt_setup::opt_setup(const opt_setup& setup) : optimiser_setup(setup)
    {
        Z0 = setup.Z0;
        er = setup.er;
        d = setup.d;
        M = setup.M;
        Zref = setup.Zref;
        freqs = setup.freqs;
        K = setup.K;
        Z_min = setup.Z_min;
        Z_max = setup.Z_max;
        Z_at_0_2 = setup.Z_at_0_2;
        Z_at_d_2 = setup.Z_at_d_2;
        Z_at_0_3 = setup.Z_at_0_3;
        Z_at_d_3 = setup.Z_at_d_3;
        R_min = setup.R_min;
        R_max = setup.R_max;
        matching_dB = setup.matching_dB;
        isolation_dB = setup.isolation_dB;
        split = setup.split;
    }

    opt_setup::opt_setup(const nlohmann::json& j) : optimiser_setup(j)
    {
        if (j.at("setup_type") != "WPD_opt")
            throw(std::logic_error("Attempted to read a setup from a different json object"));

        Z0 = j.at("Z0").get<double>();
        er = j.at("er").get<double>();
        d = j.at("d").get<double>();
        M = j.at("M").get<int>();
        Zref = j.at("Zref").get<double>();      

        freqs = j.at("freqs").get<std::vector<double>>();
        K = j.at("K").get<double>();
        Z_min = j.at("Z_min").get<double>();
        Z_max = j.at("Z_max").get<double>();
        Z_at_0_2 = j.at("Z_at_0_2").get<double>();
        Z_at_d_2 = j.at("Z_at_d_2").get<double>();
        Z_at_0_3 = j.at("Z_at_0_3").get<double>();
        Z_at_d_3 = j.at("Z_at_d_3").get<double>();
        R_min = j.at("R_min").get<double>();
        R_max = j.at("R_max").get<double>();
        matching_dB = j.at("matching_dB").get<double>();
        isolation_dB = j.at("isolation_dB").get<double>();
        split = j.at("split").get<std::vector<double>>();
    }

    nlohmann::json opt_setup::get_json() const
    {      
        return {
            { "json_type", "setup" },
            { "setup_type", "WPD_opt"},
            { "N", N },
            { "lb", lb },
            { "ub", ub },
            { "toll_bounds", toll_bounds },
            { "toll_z", toll_z },
            { "GBL_MAX", GBL_MAX },
            { "LCL_MAX", LCL_MAX },
            { "accepted_error", accepted_error },
            { "max_attempts", max_attempts },
            { "Z0", Z0 },
            { "er", er },
            { "d", d },
            { "M", M },
            { "Zref", Zref },            
            { "freqs", freqs },
            { "K", K },
            { "Z_min", Z_min },
            { "Z_max", Z_max },
            { "Z_at_0_2", Z_at_0_2 },
            { "Z_at_d_2", Z_at_d_2 },
            { "Z_at_0_3", Z_at_0_3 },
            { "Z_at_d_3", Z_at_d_3 },
            { "R_min", R_min },
            { "R_max", R_max },
            {"matching_dB", matching_dB},
            {"isolation_dB", isolation_dB},
            {"split", split},
        };
    }

    opt::opt(const opt_setup& setup) : optimiser(setup)
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
        m_Z_at_0_2 = setup.Z_at_0_2;
        m_Z_at_d_2 = setup.Z_at_d_2;
        m_Z_at_0_3 = setup.Z_at_0_3;
        m_Z_at_d_3 = setup.Z_at_d_3;
        m_R_min = setup.R_min;
        m_R_max = setup.R_max;
        m_matching_dB = setup.matching_dB;
        m_isolation_dB = setup.isolation_dB;
        m_split = setup.split;

        if (m_M > m_N / 2)
            throw(std::invalid_argument("Sine terms must be fewer or equal to the number of terms"));
        if (m_split.size() < m_freqs.size() && m_split.size() != 1)
            throw(std::invalid_argument("Number of splits & frequency points mismatch"));
        if (m_split.size() == 1)
        {
            m_split.resize(m_freqs.size());
            for (int i = 1; i < m_split.size(); i++)
                m_split[i] = m_split[0];
        }
        if (setup.N % 2 != 1)
            throw(std::invalid_argument("The number of coefficients must be odd"));
        if (m_Z0 < 1e-6 || m_er < 1e-6 || m_d < 1e-6)
            throw(std::invalid_argument("Invalid NTL physical charachteristics"));
        if (m_freqs.empty())
            throw(std::invalid_argument("Empty frequency vector"));
        if (m_Z_min < 1e-6 || m_Z_max < 1e-6 || m_Z_max < m_Z_min)
            throw(std::invalid_argument("Invalid min/max impedance(s)"));
        if (m_R_min < 1e-6 || m_R_max < 1e-6 || m_R_max < m_R_min)
            throw(std::invalid_argument("Invalid min/max resistance"));
        if (m_Z_at_0_2 < 1e-6 || m_Z_at_d_2 < 1e-6 || m_Z_at_0_3 < 1e-6 || m_Z_at_d_3 < 1e-6)
            throw(std::invalid_argument("Invalid impedance boundary conditions"));

        m_Zl.resize(m_freqs.size());
        m_split_abs.resize(m_freqs.size());


        for(int i = 0; i < m_Zl.size(); i++)
        {
            m_split_abs[i][0] = 1.0 / std::sqrt(1 + m_split[i]);     //S12
            m_split_abs[i][1] = 1.0 / std::sqrt(1 + 1 / m_split[i]); //S13
        }

        m_N2 = m_N / 2;
        m_N3 = m_N2;

        m_lb.resize(m_N);
        m_ub.resize(m_N);
        m_lb[m_N - 1] = m_R_min;
        m_ub[m_N - 1] = m_R_max;

        m_matching_abs = std::pow(10, m_matching_dB / 20);
        m_isolation_abs = std::pow(10, m_isolation_dB / 20);            

        int hw_threads = std::thread::hardware_concurrency();
        if (hw_threads == 0) hw_threads = 4;

        int m_F = m_freqs.size();

        omp_set_num_threads(std::min<int>(m_F, hw_threads));
    }

    opt::opt(const opt_setup& setup, const NTL::NTL& output2, const NTL::NTL& output3)
        : opt(setup)
    {
        m_ntl2_out = output2;
        m_ntl3_out = output3;
    }

    opt_result opt::optimise(console mode, bool output_d)
    {
        bool out = (mode == console::active);
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (m_ntl2_out.get_d() == 0) 
        {
            //Optimise the output arms
            NTL::opt_setup output_setup;
            output_setup.N = m_N2;
            output_setup.M= 0;
            output_setup.lb.assign(m_N2, -1);
            output_setup.ub.assign(m_N2, 1);
            output_setup.toll_bounds.assign(2, 1e-6);
            output_setup.toll_z.assign(m_K, 1e-6);
            output_setup.GBL_MAX = 350000 * 2;
            output_setup.LCL_MAX = 50000 * 2;
            output_setup.accepted_error = 1e-5;
            output_setup.max_attempts = 10;
            output_setup.Z0 = m_Z0;
            output_setup.er = m_er;
            output_setup.d = 80e-3;
            output_setup.Zs = { m_Zref };
            output_setup.Zl.resize(m_freqs.size());
            for (int i = 0; i < m_freqs.size(); i++)
                output_setup.Zl[i] = m_Zref * std::sqrt(m_split[i]);
            output_setup.freqs = m_freqs;
            output_setup.K = m_K;
            output_setup.Z_min = m_Z_min;
            output_setup.Z_max = m_Z_max;
            output_setup.Z_at_0 = m_Zref;
            output_setup.Z_at_d = m_Zref;            

            {
                if (out)
                    std::cout << "\nOptimising port2 output transformer...\n";
                NTL::opt output2_opt(output_setup);
                
                if(output_d)
                    m_ntl2_out = output2_opt.optimise_d(mode).ntl;
                else
                    m_ntl2_out = output2_opt.optimise(mode).ntl;
            }

            for (int i = 0; i < m_freqs.size(); i++)
                output_setup.Zl[i] = m_Zref / std::sqrt(m_split[i]);
            {
                if (out)
                    std::cout << "\nOptimising port3 output transformer...\n";
                NTL::opt output3_opt(output_setup);                                 
                if (output_d)
                    m_ntl3_out = output3_opt.optimise_d(mode).ntl;
                else
                    m_ntl3_out = output3_opt.optimise(mode).ntl;
            }
            
        }

        for(int i = 0; i < m_Zl.size(); i++)
        m_Zl[i] =
        { m_Zref,
               NTL::calculate_Zin(m_ntl2_out, m_Zref, m_freqs[i], m_K),
               NTL::calculate_Zin(m_ntl3_out, m_Zref, m_freqs[i], m_K)
        };

        if (out)
            std::cout << "\n\nStarted WPD Optimisation...\n";

        optimiser_result result = run_optimiser(mode); 

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;

        std::vector<double> Cn = result.optimised_cn;
        std::vector<double> Cn2(Cn.begin(), Cn.begin() + m_N2);
        std::vector<double> Cn3(Cn.begin() + m_N2, Cn.begin() + m_N2 + m_N3);
        double R = Cn.back();
        
        NTL::NTL ntl2(m_Z0, m_er, m_d, Cn2, m_M);
        NTL::NTL ntl3(m_Z0, m_er, m_d, Cn3, m_M);        
        WPD wpd(ntl2, ntl3, R);      
        
        if (out)
              std::cout << "*** WPD Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
  
        
        return { result, wpd, m_ntl2_out, m_ntl3_out };
    }

    void opt::print_config() const
    {
        std::cout << "========== WPD OPTIMISER CONFIGURATION ==========\n";

        std::cout << "[Physical Constraints]\n";
        std::cout << "  Z0:      " << m_Z0 << " Ohm\n";
        std::cout << "  er:      " << m_er << "\n";
        std::cout << "  d:       " << m_d << " m\n";
        std::cout << "  M:       " << m_M << " (Sine Terms)\n";
        std::cout << "  Zref:    " << m_Zref << " Ohm\n";

        std::cout << "\n[Search Space]\n";
        std::cout << "  N2 (Arm2 Coeffs): " << m_N2 << "\n";
        std::cout << "  N3 (Arm3 Coeffs): " << m_N3 << "\n";
        std::cout << "  Z Range: [" << m_Z_min << ", " << m_Z_max << "] Ohm\n";
        std::cout << "  R Range: [" << m_R_min << ", " << m_R_max << "] Ohm\n";

        std::cout << "\n[Objectives]\n";
        std::cout << "  Matching Target:  " << m_matching_dB << " dB (Linear: " << m_matching_abs << ")\n";
        std::cout << "  Isolation Target: " << m_isolation_dB << " dB (Linear: " << m_isolation_abs << ")\n";

        std::cout << "\n[Frequencies & Split]\n";
        std::cout << "  Points: " << m_freqs.size() << "\n";
        std::cout << "  K (integration steps): " << m_K << "\n";

        // Print table of Frequency vs Split Targets
        std::cout << "  Index | Freq (GHz) | Split Ratio | Target S12 (abs) | Target S13 (abs)\n";
        std::cout << "  ---------------------------------------------------------------------\n";
        for (size_t i = 0; i < m_freqs.size(); ++i)
        {
            double f_ghz = m_freqs[i] / 1e9;
            std::cout << "  " << i << "     | "
                << f_ghz << "        | "
                << m_split[i] << "         | ";

            if (i < m_split_abs.size())
                std::cout << m_split_abs[i][0] << "           | " << m_split_abs[i][1];
            else
                std::cout << "N/A";

            std::cout << "\n";
        }

        std::cout << "\n[Output Transformers (Pre-Solved)]\n";
        auto print_ntl_summary = [](const std::string& name, const NTL::NTL& ntl) {
            std::cout << "  " << name << ": ";
            if (ntl.get_d() > 0)
                std::cout << "Active (d=" << ntl.get_d() << ", N=" << ntl.get_Cn().size() << ")\n";
            else
                std::cout << "None/Passthrough\n";
            };
        print_ntl_summary("NTL2 Out", m_ntl2_out);
        print_ntl_summary("NTL3 Out", m_ntl3_out);

        std::cout << "\n[Load Impedances (Zl) - derived]\n";
        // Just print the first few to avoid screen spam if many frequencies
        int max_print = std::min<int>(5, m_Zl.size());
        for (int i = 0; i < max_print; ++i)
        {
            std::cout << "  f[" << i << "]: "
                << "Zl1=" << m_Zl[i][0] << ", "
                << "Zl2=" << m_Zl[i][1] << ", "
                << "Zl3=" << m_Zl[i][2] << "\n";
        }
        if (m_Zl.size() > max_print)
            std::cout << "  ... (" << m_Zl.size() - max_print << " more hidden)\n";

        std::cout << "=================================================\n" << std::endl;
    }

    double opt::min_objective(const std::vector<double>& Cn) const
    {
        double sum_squares{};
        

        #pragma omp parallel
        {
            double thread_sum_squares{};
            #pragma omp for nowait
            for (int i = 0; i < m_freqs.size(); i++)
            {
                double f = m_freqs[i];
                
                std::vector<double> Cn2(Cn.cbegin(), Cn.cbegin() + m_N2);
                std::vector<double> Cn3(Cn.cbegin() + m_N2, Cn.cbegin() + m_N2 + m_N3);
                double R = Cn.back();

                matrix3x3cd S = calculate_S_matrix(m_Z0, m_er, m_d, Cn2, m_d, Cn3, m_M, R, f, m_Zl[i], m_K);


                //Matching

                double S11 = std::abs(S(0, 0));
                double S22 = std::abs(S(1, 1));
                double S33 = std::abs(S(2, 2));

                if (S11 > m_matching_abs)
                    thread_sum_squares += std::pow(S11 - m_matching_abs, 2);
                if (S22 > m_matching_abs)
                    thread_sum_squares += std::pow(S22 - m_matching_abs, 2);
                if (S33 > m_matching_abs)
                    thread_sum_squares += std::pow(S33 - m_matching_abs, 2);

                //Isolation

                double S23 = std::abs(S(1, 2));
                if (S23 > m_isolation_abs)
                    thread_sum_squares += 2 * std::pow(S23 - m_isolation_abs, 2);

                //Split

                double S12 = std::abs(S(0, 1));
                double S13 = std::abs(S(0, 2));

                thread_sum_squares += 2 * std::pow(S12 - m_split_abs[i][0], 2);
                thread_sum_squares += 2 * std::pow(S13 - m_split_abs[i][1], 2);
                   
            }

            #pragma omp critical
            {
                sum_squares += thread_sum_squares;
            }
        }

        return std::sqrt(sum_squares);
    }
    void opt::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
    {
        // 1. Separate the coefficients
        const double* Cn2 = Cn;              // NTL2 starts at index 0
        const double* Cn3 = Cn + m_N2;       // NTL3 starts after NTL2

        // 2. Apply constraints (Z_start and Z_end for both arms)
        // We expect m >= 4 for a full WPD setup (Start/End for NTL2 + Start/End for NTL3)

        // --- NTL2 Constraints ---
        // Z(0) == Z0
        res[0] = NTL::calculate_Z(Cn2, m_N2, m_M, m_Z0, m_d, 0.0) - m_Z_at_0_2;
        // Z(d) == Z0
        res[1] = NTL::calculate_Z(Cn2, m_N2, m_M, m_Z0, m_d, m_d) - m_Z_at_d_3;

        // --- NTL3 Constraints ---
        // Z(0) == Z0
        res[2] = NTL::calculate_Z(Cn3, m_N3, m_M, m_Z0, m_d, 0.0) - m_Z_at_0_3;
        // Z(d) == Z0
        res[3] = NTL::calculate_Z(Cn3, m_N3, m_M, m_Z0, m_d, m_d) - m_Z_at_d_3;
    }

    void opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
    {
        // 1. Define pointers to the start of each NTL's coefficients
        const double* Cn2 = Cn;              // NTL2 starts at index 0
        const double* Cn3 = Cn + m_N2;       // NTL3 starts after NTL2 (m_N2 is size of first arm's coeffs)

        // 2. Determine points per arm
        // We assume 'm' (total constraints) is split equally between the two arms.
        unsigned k_points = m / 2;
        double dz = m_d / k_points;

        // 3. Apply constraints for NTL2 (First half of res)
        for (unsigned i = 0; i < k_points; ++i)
        {
            double z = (i + 0.5) * dz;
            // precise call: pass only m_N2 as the number of coefficients
            double Z_val = NTL::calculate_Z(Cn2, m_N2, m_M, m_Z0, m_d, z);
            res[i] = Z_val - m_Z_max;
        }

        // 4. Apply constraints for NTL3 (Second half of res)
        for (unsigned i = 0; i < k_points; ++i)
        {
            double z = (i + 0.5) * dz;
            // precise call: pass only m_N3 as the number of coefficients
            double Z_val = NTL::calculate_Z(Cn3, m_N3, m_M, m_Z0, m_d, z);
            res[k_points + i] = Z_val - m_Z_max;
        }
    }

    void opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
    {
        const double* Cn2 = Cn;
        const double* Cn3 = Cn + m_N2;

        unsigned k_points = m / 2;
        double dz = m_d / k_points;

        // NTL2
        for (unsigned i = 0; i < k_points; ++i)
        {
            double z = (i + 0.5) * dz;
            double Z_val = NTL::calculate_Z(Cn2, m_N2, m_M, m_Z0, m_d, z);
            res[i] = m_Z_min - Z_val;
        }

        // NTL3
        for (unsigned i = 0; i < k_points; ++i)
        {
            double z = (i + 0.5) * dz;
            double Z_val = NTL::calculate_Z(Cn3, m_N3, m_M, m_Z0, m_d, z);
            res[k_points + i] = m_Z_min - Z_val;
        }
    }

    // [Inside WPD_opt class implementation]

    double opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
    {
        double total_sum_squares = 0.0;
        bool calc_grad = !grad.empty();

        if (calc_grad)
            std::fill(grad.begin(), grad.end(), 0.0);

        // Unpack Inputs
        // Cn structure: [ Cn2 (m_N2) | Cn3 (m_N3) | R (1) ]
        std::vector<double> Cn2(Cn.begin(), Cn.begin() + m_N2);
        std::vector<double> Cn3(Cn.begin() + m_N2, Cn.begin() + m_N2 + m_N3);
        double R = Cn.back();

        int total_params = Cn.size();

        // Parallelize over frequencies
        #pragma omp parallel
        {
            double thread_sum_squares = 0.0;
            std::vector<double> thread_grad(total_params, 0.0);

            #pragma omp for nowait
            for (int i = 0; i < m_freqs.size(); i++)
            {
                double f = m_freqs[i];

                // 1. Calculate S and Gradient for this frequency
                auto [S, dS] = calculate_S_matrix_with_grad(
                    m_Z0, m_er, m_d, Cn2, m_d, Cn3, m_M, R, f, m_Zl[i], m_K);

                // Helper to accumulate error and gradient for a specific S-parameter target
                // Objective term: E = ( |S_xy| - Target )^2
                // dE/dCn = 2 * ( |S_xy| - Target ) * d|S_xy|/dCn
                // d|S|/dx = (1/|S|) * Re( S* . dS/dx )
                auto add_objective = [&](int r, int c, double target, double weight)
                    {
                        std::complex<double> val = S(r, c);
                        double mag = std::abs(val);

                        // Penalize only if magnitude exceeds target (for matching/isolation)
                        // or match target exactly (for split) depending on logic.
                        // Based on min_objective, Split is exact match, others are inequalities.

                        double diff = 0.0;
                        bool is_inequality = false; // Is it essentially S < Target?

                        // Identify type based on target value relative to constraints logic
                        // Using exact logic from min_objective:

                        // Matching (S11, S22, S33) -> Inequality (Minimize if > max)
                        if ((r == 0 && c == 0) || (r == 1 && c == 1) || (r == 2 && c == 2)) {
                            if (mag > m_matching_abs) {
                                diff = mag - m_matching_abs;
                                is_inequality = true;
                            }
                        }
                        // Isolation (S23) -> Inequality
                        else if (r == 1 && c == 2) {
                            if (mag > m_isolation_abs) {
                                diff = mag - m_isolation_abs;
                                is_inequality = true;
                            }
                        }
                        // Split (S12, S13) -> Equality (Match target)
                        else {
                            // Note: m_split_abs indexes map to freq index i
                            // m_split_abs[i][0] is S12, [1] is S13
                            double target_val = (c == 1) ? m_split_abs[i][0] : m_split_abs[i][1];
                            diff = mag - target_val;
                        }

                        if (std::abs(diff) > 1e-15) // Only accumulate if there is error
                        {
                            thread_sum_squares += weight * diff * diff;

                            if (calc_grad)
                            {
                                // Factor: 2 * weight * diff * (1/mag)
                                // If mag is tiny, derivative is unstable, handle gracefully
                                if (mag > 1e-12)
                                {
                                    double pre_factor = 2.0 * weight * diff / mag;
                                    for (int p = 0; p < total_params; ++p)
                                    {
                                        std::complex<double> dVal = dS[p](r, c);
                                        // Re( S* . dS/dx )
                                        double dMag = std::real(std::conj(val) * dVal);
                                        thread_grad[p] += pre_factor * dMag;
                                    }
                                }
                            }
                        }
                    };

                // Apply Objectives
                // Matching
                add_objective(0, 0, m_matching_abs, 1.0); // S11
                add_objective(1, 1, m_matching_abs, 1.0); // S22
                add_objective(2, 2, m_matching_abs, 1.0); // S33

                // Isolation (weight 2 from min_objective)
                add_objective(1, 2, m_isolation_abs, 2.0); // S23

                // Split (weight 2 from min_objective)
                add_objective(0, 1, 0.0, 2); // S12 (Target handled inside lambda)
                add_objective(0, 2, 0.0, 2); // S13 (Target handled inside lambda)
            }

            #pragma omp critical
            {
                total_sum_squares += thread_sum_squares;
                if (calc_grad)
                {
                    for (int p = 0; p < total_params; ++p)
                        grad[p] += thread_grad[p];
                }
            }
        }

        return std::sqrt(total_sum_squares);
    }
}