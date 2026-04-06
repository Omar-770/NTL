#include "flt_opt.h"
#include "common/helpers.h"

namespace flt
{
	opt_setup::opt_setup(const opt_setup& setup) : optimiser_setup(setup)
	{
		Z0 = setup.Z0;
		er = setup.er;
		d = setup.d;
		M = setup.M;
		Zs = setup.Zs;
		Zl = setup.Zl;
		freqs = setup.freqs;
		K = setup.K;
		Z_min = setup.Z_min;
		Z_max = setup.Z_max;
		Z_at_0 = setup.Z_at_0;
		Z_at_d = setup.Z_at_d;
	}

	opt_setup::opt_setup(const optimiser_setup& setup) : optimiser_setup(setup)
	{
	}

	opt_setup::opt_setup(const nlohmann::json& j) : optimiser_setup(j)
	{
		if (j.at("setup_type") != "NTL_opt")
			throw(std::logic_error("Attempted to read a setup from a different json object"));

		Z0 = j.at("Z0").get<double>();
		er = j.at("er").get<double>();
		d = j.at("d").get<double>();
		M = j.at("M").get<double>();
		Zs = j.at("Zs").get<std::vector<std::complex<double>>>();
		Zl = j.at("Zl").get<std::vector<std::complex<double>>>();
		freqs = j.at("freqs").get<std::vector<double>>();
		K = j.at("K").get<double>();
		Z_min = j.at("Z_min").get<double>();
		Z_max = j.at("Z_max").get<double>();
		Z_at_0 = j.at("Z_at_0").get<double>();
		Z_at_d = j.at("Z_at_d").get<double>();
	}

	nlohmann::json opt_setup::get_json() const
	{
		return {
			{ "json_type", "setup" },
			{ "setup_type", "NTL_opt"},
			{ "N", N },
			{ "GBL_MAX", GBL_MAX },
			{ "LCL_MAX", LCL_MAX },
			{ "accepted_error", accepted_error },
			{ "max_attempts", max_attempts },
			{ "Z0", Z0 },
			{ "er", er },
			{ "d", d },
			{ "M", M },
			{ "Zs", Zs },
			{ "Zl", Zl },
			{ "freqs", freqs },
			{ "K", K },
			{ "Z_min", Z_min },
			{ "Z_max", Z_max },
			{ "Z_at_0", Z_at_0 },
			{ "Z_at_0", Z_at_d }
		};
	}

	opt::opt(const opt_setup& setup) : optimiser(setup)
	{
		m_Z0 = setup.Z0;
		m_er = setup.er;
		m_d = setup.d;
		m_M = setup.M;
		m_Zs = setup.Zs;
		m_Zl = setup.Zl;
		m_freqs = setup.freqs;
		m_K = setup.K;
		m_Z_min = setup.Z_min;
		m_Z_max = setup.Z_max;
		m_Z_at_0 = setup.Z_at_0;
		m_Z_at_d = setup.Z_at_d;

		if (m_M > m_N)
			throw(std::invalid_argument("Sine terms must be fewer or equal to the number of terms"));

		if (m_Zl.size() != 1 && m_Zl.size() != m_freqs.size())
			throw std::invalid_argument("Zl size must be 1 or equal to frequencies count");
		if (m_Zs.size() != 1 && m_Zs.size() != m_freqs.size())
			throw std::invalid_argument("Zs size must be 1 or equal to frequencies count");
		if (m_Zl.size() == 1)
		{
			m_Zl.resize(m_freqs.size());
			for (int i = 1; i < m_Zl.size(); i++)
				m_Zl[i] = m_Zl[0];
		}
		if (m_Zs.size() == 1)
		{
			m_Zs.resize(m_freqs.size());
			for (int i = 1; i < m_Zs.size(); i++)
				m_Zs[i] = m_Zs[0];
		}

		if (m_Z0 < 1e-6 || m_er < 1e-6 || m_d < 1e-6)
			throw(std::invalid_argument("Invalid NTL physical charachteristics"));
		for (auto& z : m_Zs)
			if (std::abs(z) < 1e-6)
				throw(std::invalid_argument("Invalid source impedance(s)"));
		for (auto& z : m_Zl)
			if (std::abs(z) < 1e-6)
				throw(std::invalid_argument("Invalid load impedance(s)"));
		if (m_freqs.empty())
			throw(std::invalid_argument("Empty frequency vector"));
		if (m_Z_min < 1e-6 || m_Z_max < 1e-6 || m_Z_max < m_Z_min)
			throw(std::invalid_argument("Invalid min/max impedance(s)"));
		if (m_Z_at_0 < 1e-6 || m_Z_at_d < 1e-6)
			throw(std::invalid_argument("Invalid impedance boundary conditions"));

		if (m_lb.empty() || m_lb.size() != m_N)
			m_lb.assign(m_N, -1.0);
		if (m_ub.empty() || m_ub.size() != m_N)
			m_ub.assign(m_N, 1.0);

		if (m_toll_bounds.empty())
			m_toll_bounds.assign(2, 1e-6);
		if (m_toll_z.empty())
			m_toll_z.assign(m_K, 1e-6);

		omp_set_num_threads(std::min<int>(m_freqs.size(), 11));
	}

	opt_result opt::optimise(console mode)
	{
		return opt_result();
	}
}