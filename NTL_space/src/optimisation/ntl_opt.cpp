#include "NTL_opt.h"
#include <iostream>

namespace NTL
{
	NTL NTL_opt::optimise(bool output)
	{
		double overall_best_error = std::numeric_limits<double>::max();
		std::vector<double> overall_best_Cn(m_N);
		int attempt = 1;


		while (true)
		{
			if(output)
			{
				std::cout << "\n============================================================" << std::endl;
				std::cout << "               STARTING OPTIMIZATION ATTEMPT #" << attempt << std::endl;
				std::cout << "============================================================" << std::endl;
			}

			std::vector<double> Cn_this_attempt(m_N);
			std::vector<double> best_Cn_this_attempt(m_N);
			double global_minf, final_minf_this_attempt;

			// --- 1. GLOBAL SEARCH PHASE ---
			if(output)
				std::cout << "--- Starting Global Search Phase (using GN_ISRES) ---" << std::endl;
			nlopt::opt global_optimizer(nlopt::GN_ISRES, m_N);
			global_optimizer.set_lower_bounds(m_lb);
			global_optimizer.set_upper_bounds(m_ub);

			global_optimizer.set_min_objective([](const std::vector<double>& Cn, std::vector<double>& grad, void* data) -> double {
				return (static_cast<NTL_opt*>(data))->min_objective(Cn);
				}, this);
			global_optimizer.add_equality_mconstraint([](unsigned m, double* res, unsigned n, const double* x, double* grad, void* data) {
				(*static_cast<NTL_opt*>(data)).equality_constraints(m, res, n, x);
				}, this, m_toll_bounds);
			global_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
				(*static_cast<NTL_opt*>(data)).inequality_constraints_Zmax(m, res, n, Cn);
				}, this, m_toll_z);
			global_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
				(*static_cast<NTL_opt*>(data)).inequality_constraints_Zmin(m, res, n, Cn);
				}, this, m_toll_z);

			global_optimizer.set_ftol_rel(1e-3);
			global_optimizer.set_maxeval(m_GBL_MAX);

			try
			{
				global_optimizer.optimize(Cn_this_attempt, global_minf);
				if (output)
					std::cout << "Global search finished. Candidate error: " << global_minf << std::endl;
			}
			catch (const std::exception& ex)
			{
				std::cerr << "Global optimizer failed on attempt #" << attempt << ": " << ex.what() << std::endl;
				attempt++;
				continue;
			}

			// --- 2. LOCAL REFINEMENT PHASE ---
			if (output)
				std::cout << "\n--- Starting Local Refinement Phase (using LD_MMA) ---" << std::endl;
			nlopt::opt local_optimizer(nlopt::LD_AUGLAG, m_N);
			//nlopt::opt local_optimizer(nlopt::LN_COBYLA, N);
			local_optimizer.set_lower_bounds(m_lb);
			local_optimizer.set_upper_bounds(m_ub);


			nlopt::opt inner_optimizer(nlopt::LD_MMA, m_N);
			inner_optimizer.set_ftol_rel(1e-7);
			inner_optimizer.set_xtol_rel(1e-7);

			local_optimizer.set_local_optimizer(inner_optimizer);
			local_optimizer.set_min_objective([](const std::vector<double>& Cn, std::vector<double>& grad, void* data) -> double {
				return (static_cast<NTL_opt*>(data))->objective_with_fd_gradient(Cn, grad, data);
				}, this);
			local_optimizer.add_equality_mconstraint([](unsigned m, double* res, unsigned n, const double* x, double* grad, void* data) {
				(*static_cast<NTL_opt*>(data)).equality_constraints(m, res, n, x);
				}, this, m_toll_bounds);
			local_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
				(*static_cast<NTL_opt*>(data)).inequality_constraints_Zmax(m, res, n, Cn);
				}, this, m_toll_z);
			local_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
				(*static_cast<NTL_opt*>(data)).inequality_constraints_Zmin(m, res, n, Cn);
				}, this, m_toll_z);

			local_optimizer.set_ftol_rel(1e-7);
			local_optimizer.set_xtol_rel(1e-7);
			local_optimizer.set_maxeval(m_LCL_MAX);

			best_Cn_this_attempt = Cn_this_attempt;

			try
			{
				local_optimizer.optimize(best_Cn_this_attempt, final_minf_this_attempt);
				if (output)
					std::cout << "Local refinement finished. Error for this attempt: " << final_minf_this_attempt << std::endl;
			}
			catch (const std::exception& ex)
			{
				std::cerr << "Local optimizer failed on attempt #" << attempt << ": " << ex.what() << std::endl;
				attempt++;
				continue;
			}


			if (final_minf_this_attempt < overall_best_error)
			{
				overall_best_error = final_minf_this_attempt;
				overall_best_Cn = best_Cn_this_attempt;
				if (output)
				{
					std::cout << "\n  *************************************" << std::endl;
					std::cout << "  ** New best solution found! Error: " << overall_best_error << "\t at " << std::endl;
					std::cout << overall_best_Cn << std::endl;
					std::cout << "  *************************************" << std::endl;
				}

			}

			if (final_minf_this_attempt < m_accepted_error)
			{
				if (output)
					std::cout << "\nSUCCESS: Solution found with error " << final_minf_this_attempt
						<< ", which is below the acceptance threshold of " << m_accepted_error << "." << std::endl;
				break;
			}

			attempt++;
		}

		/// RESULTS ///
		if (output)
		{
			std::cout << "\n\n\n==============================================" << std::endl;
			std::cout << "           OPTIMIZATION COMPLETE" << std::endl;
			std::cout << "==============================================" << std::endl;
			std::cout << "Final best error found: " << overall_best_error << "\t at ";
			std::cout << overall_best_Cn << std::endl;
		}


		return NTL(m_Z0, m_er, m_d, overall_best_Cn);
	}

	/// NLOPT_FUNCTIONS

	double NTL_opt::min_objective(const std::vector<double>& Cn) const
	{
		double sum_squares{};

		for (double f : m_freqs)
		{
			matrix2x2cd T = calculate_T_matrix(m_Z0, m_er, m_d, Cn, f, m_K);

			std::complex<double> a = T(0, 0), b = T(0, 1), c = T(1, 0), d = T(1, 1);
			std::complex<double> Zin = (m_Zl * a + b) / (m_Zl * c + d);
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