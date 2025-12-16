#pragma once
#include "optimisation/optimiser.h"
#include "optimisation/WPD_opt.h"
#include "common/file_handler.h"
#include <memory>
#include "common/helpers.h"
#include <thread>

namespace WPD
{
	class opt_2;

	class opt_2 : public optimiser
	{
	public:
		opt_2(const opt_setup& setup);

		opt_result optimise(console mode = console::active);

		void print_result_logs();

	public:
		void optimise_d_arms() { m_d_arms = true; }
		void optimise_d_outputs() { m_d_outputs = true; }

	private:
		void optimise_arms();
		void optimise_transformers();
		void optimise_R();

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
		double m_R_min;
		double m_R_max;
		std::vector<double> m_split;

	private:
		//Internal common variables
		bool m_out;
		int m_F;
		int m_N1;
		bool m_d_arms;
		bool m_d_outputs;

		std::vector<std::complex<double>> m_arm2_Zl;
		std::vector<std::complex<double>> m_arm2_Zs;
		std::vector<std::complex<double>> m_arm3_Zl;
		std::vector<std::complex<double>> m_arm3_Zs;

	private:
		//Final Result
		WPD m_wpd;
		std::array<NTL::NTL, 4> m_ntl;
		double m_R;

	private:
		// nlopt functions
		double min_objective(const std::vector<double>& Cn) const override;
		void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const override;
		double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const override;
	};
}