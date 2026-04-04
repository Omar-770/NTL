#pragma once

#include "optimiser.h"
#include "models/tj.h"
#include "models/ntl.h"
#include "optimisation/NTL_opt.h"
#include <iostream>
#include <nlohmann/json.hpp>
#include <chrono>
#include <common/file_handler.h>

namespace TJ
{
	class opt;
	struct opt_setup;
	struct opt_result;

	using NTL::optimiser;
	using NTL::optimiser_setup;
	using NTL::optimiser_result;
	using NTL::console;

	struct opt_setup : public optimiser_setup
	{
		opt_setup() {};
		opt_setup(const opt_setup& setup);
		opt_setup(const nlohmann::json& j);

		double Z0{ 0 };
		double er{ 0 };
		double d{ 0 };
		int M{ 0 };
		double Zref{ 0 };
		std::vector<double> freqs{};
		int K{ 0 };
		double Z_min{ 0 };
		double Z_max{ 0 };
		std::vector<double> split{};

		nlohmann::json get_json() const override;
	};

	struct opt_result : public optimiser_result
	{
		TJ tj;
		NTL::NTL arm2;
		NTL::NTL arm3;
		NTL::opt_setup arm2_setup;
		NTL::opt_setup arm3_setup;
	};

	class opt : public optimiser
	{
	public:
		opt(const opt_setup& setup);

		opt_result optimise(console mode = console::active);

		void print_result_logs();
		void print_debug_logs();

	public:
		void optimise_d_arms() { m_opt_d_arms = true; }

	private:
		void optimise_arms();

	private:
		//General parameters
		double m_Z0;
		double m_er;
		double m_d;
		double m_Zref;
		int m_M;
		std::vector<double> m_freqs;
		int m_K;
		double m_Z_min;
		double m_Z_max;
		std::vector<double> m_split;

	private:
		//Internal common variables
		bool m_out;
		int m_F;
		bool m_opt_d_arms;

		std::vector<std::complex<double>> m_arm2_Zl;
		std::vector<std::complex<double>> m_arm2_Zs;
		std::vector<std::complex<double>> m_arm3_Zl;
		std::vector<std::complex<double>> m_arm3_Zs;

	private:
		//Final Result
		TJ m_tj;
		std::array<NTL::NTL, 2> m_ntl;
		NTL::opt_setup m_arm2_setup;
		NTL::opt_setup m_arm3_setup;

	private:
		// nlopt functions (inherited but unused since we delegate to NTL::opt)
		double min_objective(const std::vector<double>& Cn) const override;
		void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const override;
		double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const override;
	};

}