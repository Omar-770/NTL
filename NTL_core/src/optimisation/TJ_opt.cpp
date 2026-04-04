#include "TJ_opt.h"
#include "common/helpers.h"

namespace TJ
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
		split = setup.split;
	}

	opt_setup::opt_setup(const nlohmann::json& j) : optimiser_setup(j)
	{
		if (j.at("setup_type") != "TJ_opt")
			throw(std::logic_error("Attempted to read a setup from a different json object"));

		Z0 = j.at("Z0").get<double>();
		er = j.at("er").get<double>();
		d = j.at("d").get<double>();
		M = j.at("M").get<int>();
		Zref = j.at("Zref").get<double>();
		freqs = j.at("freqs").get<std::vector<double>>();
		K = j.at("K").get<int>();
		Z_min = j.at("Z_min").get<double>();
		Z_max = j.at("Z_max").get<double>();
		split = j.at("split").get<std::vector<double>>();
	}

	nlohmann::json opt_setup::get_json() const
	{
		return {
			{ "json_type", "setup" },
			{ "setup_type", "TJ_opt"},
			{ "N", N },
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
			throw(std::invalid_argument("Invalid NTL physical characteristics"));
		if (m_freqs.empty())
			throw(std::invalid_argument("Empty frequency vector"));
		if (m_Z_min < 1e-6 || m_Z_max < 1e-6 || m_Z_max < m_Z_min)
			throw(std::invalid_argument("Invalid min/max impedance(s)"));

		m_F = m_freqs.size();

		auto set_physical_params = [&](NTL::NTL& ntl)
			{
				ntl.set_Z0(m_Z0);
				ntl.set_er(m_er);
				ntl.set_d(m_d);
				ntl.set_M(m_M);
			};

		for (auto& ntl : m_ntl)
			set_physical_params(ntl);

		m_tj.set_Z0(m_Z0);
		m_tj.set_er(m_er);

		auto resize_vec = [&](std::vector<std::complex<double>>& vec)
			{
				vec.resize(m_F);
			};

		resize_vec(m_arm2_Zl);
		resize_vec(m_arm2_Zs);
		resize_vec(m_arm3_Zl);
		resize_vec(m_arm3_Zs);

		m_opt_d_arms = false;
	}

	opt_result opt::optimise(console mode)
	{
		m_out = mode == console::active ? true : false;
		m_use_local_grad_free = true;

		if (m_out)
		{
			std::cout << "Initialising T-Junction optimisation...\n";
			std::cout << "Optimising arms independently using NTL::opt...\n\n";
		}

		optimise_arms();

		m_tj.set_arms(m_ntl[0], m_ntl[1]);

		return opt_result{ {}, m_tj, m_ntl[0], m_ntl[1], m_arm2_setup, m_arm3_setup };
	}

	void opt::print_result_logs()
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
			std::cout << "NTL3: (Gamma) = " << m_ntl[1].S11(f, m_arm3_Zs[i], m_arm3_Zl[i]) << "\n\n";
		}

		std::cout << "m_arm2_Zs: " << m_arm2_Zs << std::endl;
		std::cout << "m_arm2_Zl: " << m_arm2_Zl << std::endl;
		std::cout << "m_arm3_Zs: " << m_arm3_Zs << std::endl;
		std::cout << "m_arm3_Zl: " << m_arm3_Zl << std::endl;

		std::cout << "Final TJ:\n";
		std::cout << "Arm2: (Cn) = " << m_ntl[0].get_Cn() << '\n';
		std::cout << "Arm3: (Cn) = " << m_ntl[1].get_Cn() << '\n';
	}

	void opt::print_debug_logs()
	{
		std::cout << "\n========================================\n";
		std::cout << "       TJ_OPT DEBUG LOGS                \n";
		std::cout << "========================================\n";

		std::cout << "-- Physical Constants --\n";
		std::cout << "Z0: " << m_Z0 << " Ohm\n";
		std::cout << "Epsilon_r: " << m_er << "\n";
		std::cout << "Fourier Sine Terms (M): " << m_M << "\n";
		std::cout << "Integration Steps (K): " << m_K << "\n";

		std::cout << "\n-- Current Global Variables --\n";
		std::cout << "Arm Length (m_d): " << m_d << " m (" << m_d * 1000.0 << " mm)\n";
		std::cout << "Reference Z (m_Zref): " << m_Zref << " Ohm\n";

		std::cout << "\n-- Optimization Constraints --\n";
		std::cout << "Z_min: " << m_Z_min << " Ohm\n";
		std::cout << "Z_max: " << m_Z_max << " Ohm\n";
		std::cout << "Current Coeffs Count (m_N): " << m_N << "\n";

		std::cout << "\n-- Configured Frequencies & Splits --\n";
		std::cout << "Count (m_F): " << m_F << "\n";
		for (size_t i = 0; i < m_freqs.size(); ++i)
		{
			std::cout << "  Freq " << i << ": " << m_freqs[i] / 1e9
				<< " GHz \t| Split: " << m_split[i] << "\n";
		}

		std::cout << "\n-- Calculated Target Impedances (Per Frequency) --\n";
		std::cout << "These are the impedances the arms must match:\n";
		for (size_t i = 0; i < m_F; ++i)
		{
			std::cout << "  [" << m_freqs[i] / 1e9 << " GHz] "
				<< "Arm2_Zs: " << m_arm2_Zs[i] << " | "
				<< "Arm2_Zl: " << m_arm2_Zl[i] << "\n";
			std::cout << "             "
				<< "Arm3_Zs: " << m_arm3_Zs[i] << " | "
				<< "Arm3_Zl: " << m_arm3_Zl[i] << "\n";
		}

		std::cout << "\n-- Internal NTL Component States --\n";

		auto print_ntl = [](const char* name, NTL::NTL& n) {
			std::cout << "[" << name << "]\n";
			std::cout << "  Length: " << n.get_d() * 1000.0 << " mm\n";
			std::cout << "  Coeffs (Cn): " << n.get_Cn() << "\n";
			std::cout << "  Start Z(0): " << n.Z(0.0) << " Ohm\n";
			std::cout << "  End Z(d):   " << n.Z(n.get_d()) << " Ohm\n";
			};

		print_ntl("NTL 0 (Arm 2)", m_ntl[0]);
		print_ntl("NTL 1 (Arm 3)", m_ntl[1]);

		std::cout << "\n-- Flags --\n";
		std::cout << "Optimize Arm Length (m_d_arms): " << (m_opt_d_arms ? "TRUE" : "FALSE") << "\n";

		std::cout << "========================================\n\n";
	}

	void opt::optimise_arms()
	{
		NTL::opt_setup setup;

		setup.N = m_N;
		setup.lb.assign(m_N, -1);
		setup.ub.assign(m_N, 1);
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

		setup.Zs.resize(m_F);
		setup.Zl.resize(m_F);

		// ================= Arm 2 =================
		if (m_out) std::cout << "\n\t-->Optimising first arm (Arm 2)...\n\n";

		for (int i = 0; i < m_F; i++) {
			m_arm2_Zs[i] = m_Zref * (1.0 + m_split[i]);
			m_arm2_Zl[i] = m_Zref;

			setup.Zs[i] = m_arm2_Zs[i];
			setup.Zl[i] = m_arm2_Zl[i];		}



		NTL::opt opt2(setup);
		m_arm2_setup = setup;
		if (m_opt_d_arms)
			m_ntl[0] = opt2.optimise_d(m_out ? console::active : console::inactive).ntl;
		else
			m_ntl[0] = opt2.optimise(m_out ? console::active : console::inactive).ntl;


		// ================= Arm 3 =================
		if (m_out) std::cout << "\n\t-->Optimising second arm (Arm 3)...\n\n";

		for (int i = 0; i < m_F; i++) {
			m_arm3_Zs[i] = m_Zref * (1.0 + 1.0 / m_split[i]);
			m_arm3_Zl[i] = m_Zref;

			setup.Zs[i] = m_arm3_Zs[i];
			setup.Zl[i] = m_arm3_Zl[i];
		}

		setup.Z_at_0 = m_Zref * (1.0 + 1.0 / m_split[0]);
		setup.Z_at_d = m_Zref;

		NTL::opt opt3(setup);
		m_arm3_setup = setup;
		if (m_opt_d_arms)
			m_ntl[1] = opt3.optimise_d(m_out ? console::active : console::inactive).ntl;
		else
			m_ntl[1] = opt3.optimise(m_out ? console::active : console::inactive).ntl;
	}

	// ================= Dummy overrides =================
	// These are satisfied for the optimiser interface, but unused because 
	// the T-junction delegates directly to independent NTL::opt instances.

	double opt::min_objective(const std::vector<double>& Cn) const
	{
		return 0.0;
	}

	void opt::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
	{}

	void opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
	{}

	void opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
	{}

	double opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
	{
		return 0.0;
	}
}