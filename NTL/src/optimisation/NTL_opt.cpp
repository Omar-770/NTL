#include "NTL_opt.h"

namespace NTL
{
	opt_setup::opt_setup(const opt_setup& setup): optimiser_setup(setup)
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

		omp_set_num_threads(std::min<int>(m_freqs.size(), 11));
	}

	opt_result opt::optimise(console mode)
	{
		bool out = (mode == console::active) ? true : false;
		auto start_time = std::chrono::high_resolution_clock::now();
		optimiser_result result = run_optimiser(mode);
		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;
		NTL ntl(m_Z0, m_er, m_d, result.optimised_cn, m_M);
		

		if (out)
			std::cout << "*** Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		return opt_result(result, ntl);
	}

	opt_result opt::optimise_d(double resolution, console mode)
	{
		bool out = (mode == console::active) ? true : false;		

		double init_d = m_d;
		int init_max_attempts = m_max_attempts;

		auto start_time = std::chrono::high_resolution_clock::now();
		if (out)
			std::cout << "\n\t===>>> Starting with d = " << m_d * 1000 << "mm" << std::endl;
		optimiser_result init_result = run_optimiser(mode);
		optimiser_result new_result = init_result;

		if (init_result.final_error < m_accepted_error)
		{
			double error_this_attempt = init_result.final_error;

			while (true)
			{
				m_d -= resolution;

				if (out)
					std::cout << "\n\t===>>> Trimming NTL by " << resolution * 1000 << "mm to:\t"
					<< m_d * 1000 << "mm" << std::endl;

				optimiser_result result_this_attempt = run_optimiser(mode);

				if (result_this_attempt.final_error > m_accepted_error)
				{
					m_d += resolution;
					break;
				}

				new_result = result_this_attempt;
			}

		}

		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;

		NTL ntl(m_Z0, m_er, m_d, new_result.optimised_cn, m_M);

		if (out)
		{
			std::cout << "Final NTL length:\t" << ntl.get_d() * 1000 << "mm" << std::endl;
			std::cout << "*** Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		}

		m_d = init_d;
		m_max_attempts = init_max_attempts;

		return opt_result(new_result, ntl);
	}

	opt_result opt::optimise_d(console mode)
	{
		return optimise_d(1.0e-3, mode);
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
				std::complex<double> Zl = m_Zl[i];
				std::complex<double> Zs = m_Zs[i];

				matrix2x2cd T = calculate_T_matrix(m_Z0, m_er, m_d, Cn, m_M, f, m_K);

				std::complex<double> a = T(0, 0), b = T(0, 1), c = T(1, 0), d = T(1, 1);
				std::complex<double> Zin = (Zl * a + b) / (Zl * c + d);
				std::complex<double> gamma = (Zin - std::conj(Zs)) / (Zin + Zs);
				thread_sum_squares += std::pow(std::norm(gamma), 2);
			}

			#pragma omp critical
			{
				sum_squares += thread_sum_squares;
			}
		}

		return sum_squares;
	}

	void opt::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		res[0] = calculate_Z(Cn, n, m_M, m_Z0, m_d, 0) - m_Z_at_0;
		res[1] = calculate_Z(Cn, n, m_M, m_Z0, m_d, m_d) - m_Z_at_d;
	}

	void opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = calculate_Z(Cn, n, m_M, m_Z0, m_d, (i + 0.5) * dz) - m_Z_max;
	}

	void opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = m_Z_min - calculate_Z(Cn, n, m_M, m_Z0, m_d, (i + 0.5) * dz);
	}

	double opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
	{
		double total_sum_squares{};
		int F = m_freqs.size();
		int N = Cn.size();

		if (!grad.empty())
			std::fill(grad.begin(), grad.end(), 0.0);
		
		#pragma omp parallel
		{
			double thread_sum_squares{};
			std::vector<double> thread_grad(N, 0);

			#pragma omp for nowait
			for (int i = 0; i < F; i++)
			{
				double f = m_freqs[i];
				std::complex<double> Zl = m_Zl[i];
				std::complex<double> Zs = m_Zs[i];

				auto [T, dT] = calculate_T_matrix_with_grad(m_Z0, m_er, m_d, Cn, m_M, f, m_K);

				std::complex<double> A = T(0, 0);
				std::complex<double> B = T(0, 1);
				std::complex<double> C = T(1, 0);
				std::complex<double> D = T(1, 1);

				std::complex<double> num = A * Zl + B;
				std::complex<double> den = C * Zl + D;

				std::complex<double> Zin = num / den;
				std::complex<double> gamma = (Zin - std::conj(Zs)) / (Zin + Zs);
				std::complex<double> multiplication_factor = 8.0 * std::norm(gamma) * Zs.real() / (den * std::pow(Zin + Zs, 2));
				thread_sum_squares += std::pow(std::norm(gamma), 2);
				

				for (int n = 0; n < N; n++)
				{
					std::complex<double> dA = dT[n](0, 0);
					std::complex<double> dB = dT[n](0, 1);
					std::complex<double> dC = dT[n](1, 0);
					std::complex<double> dD = dT[n](1, 1);
					std::complex<double> dgamma = Zl * dA + dB - Zin * (Zl * dC + dD);
					dgamma *= multiplication_factor;

					thread_grad[n] +=  (gamma.real() * dgamma.real() + gamma.imag() * dgamma.imag());

				}
			}

			#pragma omp critical
			{
				total_sum_squares += thread_sum_squares;
				if(!grad.empty())
					for (int n = 0; n < N; n++)
						grad[n] += thread_grad[n];
			}

		}

		return total_sum_squares;
	}
}