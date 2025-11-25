#include "optimiser.h"
#include <iostream>

namespace NTL
{
	opt_setup::opt_setup(const nlohmann::json& j)
	{
		N = j.at("N").get<int>(); lb = j.at("lb").get<std::vector<double>>();
		lb = j.at("lb").get<std::vector<double>>(); ub = j.at("ub").get<std::vector<double>>(); 
		toll_bounds = j.at("toll_bounds").get<std::vector<double>>();
		toll_z = j.at("toll_z").get<std::vector<double>>(); 
		GBL_MAX = j.at("GBL_MAX").get<double>();
		LCL_MAX = j.at("LCL_MAX").get<double>();
		accepted_error = j.at("accepted_error").get<double>();
		max_attempts = j.at("max_attempts").get<double>();			
	}

	nlohmann::json opt_setup::get_json() const
	{
		return {
			{ "json_type", "setup" },
			{ "setup_type", "opt"},
			{ "N", N },
			{ "lb", lb },
			{ "ub", ub },
			{ "toll_bounds", toll_bounds },
			{ "toll_z", toll_z },
			{ "GBL_MAX", GBL_MAX },
			{ "LCL_MAX", LCL_MAX },
			{ "accepted_error", accepted_error },
			{ "max_attempts", max_attempts },
		};
	}

	opt::opt(const opt_setup& setup) : m_N(setup.N), m_lb(setup.lb), m_ub(setup.ub), m_toll_bounds(setup.toll_bounds),
		m_toll_z(setup.toll_z), m_GBL_MAX(setup.GBL_MAX), m_LCL_MAX(setup.LCL_MAX), m_accepted_error(setup.accepted_error),
		m_max_attempts(setup.max_attempts)
	{
		if (m_N == 0 || m_lb.empty() || m_ub.empty() || m_toll_bounds.empty() || m_toll_z.empty())
			throw(std::invalid_argument("Incomplete optimisation setup"));
	}

	opt_result opt::optimiser(console mode)
	{
		bool output = (mode == console::active) ? true : false;

		double overall_best_error = std::numeric_limits<double>::max();
		std::vector<double> overall_best_Cn(m_N);
		int attempt = 1;


		while (attempt <= m_max_attempts)
		{
			if (output)
			{
				std::cout << "\n============================================================" << std::endl;
				std::cout << "               STARTING OPTIMISATION ATTEMPT #" << attempt << std::endl;
				std::cout << "============================================================" << std::endl;
			}

			std::vector<double> Cn_this_attempt(m_N);
			std::vector<double> best_Cn_this_attempt(m_N);
			double global_minf, final_minf_this_attempt;

			// --- 1. GLOBAL SEARCH PHASE ---
			if (output)
				std::cout << "--- Starting Global Search Phase (using GN_ISRES) ---" << std::endl;
			nlopt::opt global_optimizer(nlopt::GN_ISRES, m_N);
			global_optimizer.set_lower_bounds(m_lb);
			global_optimizer.set_upper_bounds(m_ub);

			global_optimizer.set_min_objective([](const std::vector<double>& Cn, std::vector<double>& grad, void* data) -> double {
				return (static_cast<opt*>(data))->min_objective(Cn);
				}, this);
			global_optimizer.add_equality_mconstraint([](unsigned m, double* res, unsigned n, const double* x, double* grad, void* data) {
				(*static_cast<opt*>(data)).equality_constraints(m, res, n, x);
				}, this, m_toll_bounds);
			global_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
				(*static_cast<opt*>(data)).inequality_constraints_Zmax(m, res, n, Cn);
				}, this, m_toll_z);
			global_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
				(*static_cast<opt*>(data)).inequality_constraints_Zmin(m, res, n, Cn);
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

			if (global_minf > m_accepted_error)
			{
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
					return (static_cast<opt*>(data))->objective_with_fd_gradient(Cn, grad, data);
					}, this);
				local_optimizer.add_equality_mconstraint([](unsigned m, double* res, unsigned n, const double* x, double* grad, void* data) {
					(*static_cast<opt*>(data)).equality_constraints(m, res, n, x);
					}, this, m_toll_bounds);
				local_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
					(*static_cast<opt*>(data)).inequality_constraints_Zmax(m, res, n, Cn);
					}, this, m_toll_z);
				local_optimizer.add_inequality_mconstraint([](unsigned m, double* res, unsigned n, const double* Cn, double*, void* data) {
					(*static_cast<opt*>(data)).inequality_constraints_Zmin(m, res, n, Cn);
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

			}
			else
			{
				if (output)
					std::cout << "Skipping Local Refinement Phase..." << std::endl;
				best_Cn_this_attempt = Cn_this_attempt;
				final_minf_this_attempt = global_minf;
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
			std::cout << "==============================================" << std::endl;
			std::cout << "           OPTIMIZATION COMPLETE" << std::endl;
			std::cout << "==============================================" << std::endl;
			std::cout << "Final best error found: " << overall_best_error << "\t at ";
			std::cout << overall_best_Cn << std::endl;
		}


		return  { overall_best_Cn, overall_best_error, attempt };
	}

}