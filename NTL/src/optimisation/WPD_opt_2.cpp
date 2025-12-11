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
        
            
      
        // NTL optimiser setup      
        
        m_ntl_setup.accepted_error = m_accepted_error;
        m_ntl_setup.max_attempts = m_max_attempts;
        m_ntl_setup.N = m_N;
        m_ntl_setup.GBL_MAX = m_GBL_MAX;
        m_ntl_setup.LCL_MAX = m_LCL_MAX;
        m_ntl_setup.lb = m_lb;
        m_ntl_setup.ub = m_ub;
        m_ntl_setup.toll_bounds = m_toll_bounds;
        m_ntl_setup.toll_z = m_toll_z;

        m_ntl_setup.Z0 = m_Z0;
        m_ntl_setup.d = m_d;
        m_ntl_setup.er = m_er;   
        m_ntl_setup.M = 0;
        m_ntl_setup.freqs = m_freqs;
        m_ntl_setup.K = m_K;
        m_ntl_setup.Z_min = m_Z_min;
        m_ntl_setup.Z_max = m_Z_max;
        m_ntl_setup.Z_at_0 = m_Zref;
        m_ntl_setup.Z_at_d = m_Zref;

        m_lb.assign(m_N, -1);
        m_ub.assign(m_N, 1);
        m_toll_bounds.assign(2, 1e-6);
        m_toll_z.assign(m_K, 1e-6);

        m_F = m_freqs.size();
	}

    opt_result opt_2::optimise(console mode)
    {
        /// Output transformers

        //out2 Z0 -> sqrt(split) * Zref
        //out3 Z0 -> Zref / sqrt(split)

        std::vector<std::complex<double>> out2_Zl(m_F);
        std::vector<std::complex<double>> out3_Zl(m_F);
            
        for (int i = 0; i < m_F; i++)
        {
            out2_Zl[i] = m_Zref * std::sqrt(m_split[i]);
            out3_Zl[i] = m_Zref / std::sqrt(m_split[i]);
        }

        m_ntl_setup.Zl = out2_Zl;
        m_ntl_opt[0] = std::make_unique<NTL::opt>(m_ntl_setup);
        m_out2 = m_ntl_opt[0]->optimise(mode).ntl;

        m_ntl_setup.Zl = out3_Zl;
        m_ntl_opt[1] = std::make_unique<NTL::opt>(m_ntl_setup);
        m_out3 = m_ntl_opt[1]->optimise(mode).ntl;

        /// Intermediate impedance
        std::vector<std::complex<double>> out2_Zin(m_F);
        std::vector<std::complex<double>> out3_Zin(m_F);

        for (int i = 0; i < m_F; i++)
        {
            out2_Zin[i] = std::abs(m_out2.Zin(m_Zref, m_freqs[i], m_K));
            out3_Zin[i] = std::abs(m_out3.Zin(m_Zref, m_freqs[i], m_K));
        }

        /// Arms 

        m_ntl_setup.Zl = out2_Zin;
        m_ntl_opt[2] = std::make_unique<NTL::opt>(m_ntl_setup);
        m_ntl2 = m_ntl_opt[2]->optimise(mode).ntl;

        m_ntl_setup.Zl = out3_Zin;
        m_ntl_opt[3] = std::make_unique<NTL::opt>(m_ntl_setup);
        m_ntl3 = m_ntl_opt[3]->optimise(mode).ntl;

        return opt_result();
    }
}