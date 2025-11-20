#include "NTL_opt.h"

namespace NTL
{
	NTL_opt::NTL_opt(const NTL_opt_setup& setup) : opt(setup),
		m_Z0(setup.Z0), m_er(setup.er), m_d(setup.d), m_Zs(setup.Zs), m_Zl(setup.Zl),
		m_freqs(setup.freqs), m_K(setup.K), m_Z_min(setup.m_Z_min), m_Z_max(setup.m_Z_max)
	{
		if (m_Zl.size() < m_freqs.size() && m_Zl.size() != 1)
			throw(std::invalid_argument("Number of load impedances & frequency points mismatch"));

		if (m_Zl.size() == 1)
		{
			m_Zl.resize(m_freqs.size());
			for (int i = 1; i < m_Zl.size(); i++)
				m_Zl[i] = m_Zl[0];
		}

		if (m_Z0 < 1e-6 || m_er < 1e-6 || m_d < 1e-6)
			throw(std::invalid_argument("Invalid NTL physical charachteristics"));
		if (m_Zs < 1e-6)
			throw(std::invalid_argument("Invalid source impedance"));
		for (auto& z : m_Zl)
			if (z < 1e-6)
				throw(std::invalid_argument("Invalid load impedance(s)"));

		if (m_freqs.empty())
			throw(std::invalid_argument("Empty frequency vector"));
		if (m_Z_min < 1e-6 || m_Z_max < 1e-6 || m_Z_max < m_Z_min)
			throw(std::invalid_argument("Invalid min/max impedance(s)"));
	}

	NTL_opt_result NTL_opt::optimise(console mode)
	{
		bool out = (mode == console::active) ? true : false;
		NTL ntl(m_Z0, m_er, m_d);
		auto start_time = std::chrono::high_resolution_clock::now();
		opt_result opt_result = optimiser(mode);
		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;
		ntl.set_Cn(opt_result.optimised_cn);

		if (out)
			std::cout << "*** Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		return NTL_opt_result(opt_result, ntl);
	}

	NTL_opt_result NTL_opt::optimise_d(double resolution, console mode)
	{
		bool out = (mode == console::active) ? true : false;

		NTL ntl(m_Z0, m_er, m_d);

		double init_d = m_d;
		int init_max_attempts = m_max_attempts;

		auto start_time = std::chrono::high_resolution_clock::now();
		opt_result init_result = optimiser(mode);
		opt_result new_result = init_result;

		if (init_result.final_error < m_accepted_error)
		{
			double error_this_attempt = init_result.final_error;

			while (true)
			{
				m_d -= resolution;

				if (out)
					std::cout << "\n\t===>>> Trimming NTL by " << resolution * 1000 << "mm to:\t"
					 << m_d * 1000 << "mm" << std::endl;

				opt_result result_this_attempt = optimiser(mode);

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

		ntl.set_d(m_d);
		ntl.set_Cn(new_result.optimised_cn);

		if (out)
		{
			std::cout << "Final NTL length:\t" << ntl.get_d() * 1000 << "mm" << std::endl;
			std::cout << "*** Optimisation finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		}

		m_d = init_d;
		m_max_attempts = init_max_attempts;

		return NTL_opt_result(new_result, ntl);
	}

	NTL_opt_result NTL_opt::optimise_d(console mode)
	{
		return optimise_d(1.0e-3, mode);
	}


	double NTL_opt::min_objective(const std::vector<double>& Cn) const
	{

		double sum_squares{};

		for (int i = 0; i < m_freqs.size(); i++)
		{
			double f = m_freqs[i];
			double Zl = m_Zl[i];

			matrix2x2cd T = calculate_T_matrix(m_Z0, m_er, m_d, Cn, f, m_K);

			std::complex<double> a = T(0, 0), b = T(0, 1), c = T(1, 0), d = T(1, 1);
			std::complex<double> Zin = (Zl * a + b) / (Zl * c + d);
			std::complex<double> gamma = (Zin - m_Zs) / (Zin + m_Zs);
			sum_squares += std::pow(std::norm(gamma), 2);
		}

		return std::sqrt(sum_squares);
	}

	void NTL_opt::equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		res[0] = calculate_Z(Cn, n, m_Z0, m_d, 0) - m_Z0;
		res[1] = calculate_Z(Cn, n, m_Z0, m_d, m_d) - m_Z0;
	}

	void NTL_opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = calculate_Z(Cn, n, m_Z0, m_d, (i + 0.5) * dz) - m_Z_max;
	}

	void NTL_opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = m_Z_min - calculate_Z(Cn, n, m_Z0, m_d, (i + 0.5) * dz);
	}

	double NTL_opt::objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const
	{

		const double f_base = min_objective(Cn);

		if (!grad.empty())
		{
			const double h = 1e-8;
			int N = Cn.size();

#pragma omp parallel for
			for (int i = 0; i < N; ++i)
			{
				std::vector<double> Cn_nudged = Cn;
				Cn_nudged[i] += h;

				double f_nudged = min_objective(Cn_nudged);

				grad[i] = (f_nudged - f_base) / h;
			}
		}

		return f_base;
	}
}