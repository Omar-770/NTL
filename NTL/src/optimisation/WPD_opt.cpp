#include "WPD_opt.h"

namespace WPD
{
    WPD_opt_setup::WPD_opt_setup(const WPD_opt_setup& setup) : opt_setup(setup)
    {
        Z0 = setup.Z0;
        er = setup.er;
        d = setup.d;
        Zref = setup.Zref;
        ntl2 = setup.ntl2;
        ntl3 = setup.ntl3;
        freqs = setup.freqs;
        K = setup.K;
        Z_min = setup.Z_min;
        Z_max = setup.Z_max;
        R_min = setup.R_min;
        R_max = setup.R_max;
        matching_dB = setup.matching_dB;
        isolation_dB = setup.isolation_dB;
        split = setup.split;
    }

    WPD_opt_setup::WPD_opt_setup(const nlohmann::json& j) : opt_setup(j)
    {
        if (j.at("setup_type") != "WPD_opt")
            throw(std::logic_error("Attempted to read a setup from a different json object"));

        Z0 = j.at("Z0").get<double>();
        er = j.at("er").get<double>();
        d = j.at("d").get<double>();
        Zref = j.at("Zref").get<double>();
        
        ntl2.set_Z0(j.at("ntl2-Z0").get<double>());
        ntl2.set_er(j.at("ntl2-er").get<double>());
        ntl2.set_d(j.at("ntl2-d").get<double>());
        ntl2.set_Cn(j.at("ntl2-Cn").get<std::vector<double>>());

        ntl3.set_Z0(j.at("ntl3-Z0").get<double>());
        ntl3.set_er(j.at("ntl3-er").get<double>());
        ntl3.set_d(j.at("ntl3-d").get<double>());
        ntl3.set_Cn(j.at("ntl3-Cn").get<std::vector<double>>());

        freqs = j.at("freqs").get<std::vector<double>>();
        K = j.at("K").get<double>();
        Z_min = j.at("Z_min").get<double>();
        Z_max = j.at("Z_max").get<double>();
        R_min = j.at("R_min").get<double>();
        R_max = j.at("R_max").get<double>();
        matching_dB = j.at("matching_dB").get<double>();
        isolation_dB = j.at("isolation_dB").get<double>();
        split = j.at("split").get<std::vector<double>>();
    }

    nlohmann::json WPD_opt_setup::get_json() const
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
            { "Zref", Zref },

            { "ntl2-Z0", ntl2.get_Z0()},
            { "ntl2-er", ntl2.get_er()},
            { "ntl2-d", ntl2.get_d()},
            { "ntl2-Cn", ntl2.get_Cn()},

            { "ntl3-Z0", ntl3.get_Z0()},
            { "ntl3-er", ntl3.get_er()},
            { "ntl3-d", ntl3.get_d()},
            { "ntl3-Cn", ntl3.get_Cn()},
            
            { "freqs", freqs },
            { "K", K },
            { "Z_min", Z_min },
            { "Z_max", Z_max },
            { "R_min", R_min },
            { "R_max", R_max },
            {"matching_dB", matching_dB},
            {"isolation_dB", isolation_dB},
            {"split", split},
        };
    }

    WPD_opt::WPD_opt(const WPD_opt_setup& setup) : opt(setup)
    {
        m_Z0 = setup.Z0;
        m_er = setup.er;
        m_d = setup.d;
        m_Zref = setup.Zref;
        m_ntl2 = setup.ntl2;
        m_ntl3 = setup.ntl3;
        m_freqs = setup.freqs;
        m_K = setup.K;
        m_Z_min = setup.Z_min;
        m_Z_max = setup.Z_max;
        m_R_min = setup.R_min;
        m_R_max = setup.R_max;
        m_matching_dB = setup.matching_dB;
        m_isolation_dB = setup.isolation_dB;
        m_split = setup.split;

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

        m_Zl.resize(m_freqs.size());
        m_split_abs.resize(m_freqs.size());

        for(int i = 0; i < m_Zl.size(); i++)
        {
            m_Zl[i] =
            { m_Zref,
                   NTL::calculate_Zin(m_ntl2, m_Zref, m_freqs[i], m_K),
                   NTL::calculate_Zin(m_ntl3, m_Zref, m_freqs[i], m_K)
            };

            m_split_abs[i][0] = std::sqrt(1 + m_split[i]);     //S12
            m_split_abs[i][1] = std::sqrt(1 + 1 / m_split[i]); //S13
        }

        m_N2 = m_N / 2;
        m_N3 = m_N2;

        m_lb[m_N - 1] = m_R_min;
        m_ub[m_N - 1] = m_R_max;

        m_matching_abs = std::pow(10, m_matching_dB / 20);
        m_isolation_abs = std::pow(10, m_isolation_dB / 20);            

        omp_set_num_threads(std::min<int>(m_freqs.size(), 11));
    }

    double WPD_opt::min_objective(const std::vector<double>& Cn) const
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

                matrix3x3cd S = calculate_S_matrix(m_Z0, m_er, m_d, Cn2, m_d, Cn3, R, f, m_Zl[i], m_K);


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
    void WPD_opt::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
    {
    }

    void WPD_opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
    {
    }

    void WPD_opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
    {
    }

    double WPD_opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
    {
        return 0.0;
    }
}