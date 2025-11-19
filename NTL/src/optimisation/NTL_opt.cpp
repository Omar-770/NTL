#include "NTL_opt.h"

namespace NTL
{
	NTL_opt_result NTL_opt::optimise(console mode)
	{
		if (m_Z0 == 0 || m_er == 0 || m_d == 0)
			throw(std::invalid_argument("Invalid NTL physical charachteristics"));
		if (0 < m_Zs < 1e-6)
			throw(std::invalid_argument("Invalid source impedance"));
		for (auto& z : m_Zl)
			if (0 < z < 1e-6)
				throw(std::invalid_argument("Invalid load impedance(s)"));

		if (m_freqs.empty())
			throw(std::invalid_argument("Empty frequency vector"));
		if (m_Z_min == 0 || m_Z_max == 0 || m_Z_max < m_Z_min)
			throw(std::invalid_argument("Invalid min/max impedance(s)"));


		NTL ntl(m_Z0, m_er, m_d);
		opt_result opt_result = optimiser(mode);
		ntl.set_Cn(opt_result.optimised_cn);
		
		return NTL_opt_result(opt_result, ntl);
	}

	NTL_opt_result NTL_opt::optimise_d(console mode)
	{
		bool out = (mode == console::active) ? true : false;

		NTL ntl(m_Z0, m_er, m_d);
		
		double init_d = m_d;
		int init_max_attempts = m_max_attempts;

		double init_error = optimiser(mode).final_error;
		opt_result new_result;
		
		if (init_error < m_accepted_error)
		{
			double error_this_attempt = init_error;
			
			while (true)
			{
				m_d -= 1e-3;
				//m_max_attempts = 5;

				if (out)
					std::cout << "\n\t===>>> Trimming NTL by 1mm to: \t" << m_d << std::endl;

				opt_result result_this_attempt = optimiser(mode);

				if (result_this_attempt.final_error > m_accepted_error)
				{
					m_d += 1e-3;
					break;
				}

				new_result = result_this_attempt;				
			}
			
		}

		ntl.set_d(m_d);
		ntl.set_Cn(new_result.optimised_cn);

		if (out)
			std::cout << "Final NTL length:\t" << ntl.get_d() << std::endl;

		m_d = init_d;
		m_max_attempts = init_max_attempts;

		return NTL_opt_result(new_result, ntl);
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
		res[0] = calculate_Z(Cn, n, 0.0, m_Z0, m_d) - m_Z0;
		res[1] = calculate_Z(Cn, n, m_d, m_Z0, m_d) - m_Z0;
	}

	void NTL_opt::inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = calculate_Z(Cn, n, (i + 0.5) * dz, m_Z0, m_d) - m_Z_max;
	}

	void NTL_opt::inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const
	{
		double dz = m_d / m_K;
		for (int i = 0; i < m; ++i)
			res[i] = m_Z_min - calculate_Z(Cn, n, (i + 0.5) * dz, m_Z0, m_d);
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