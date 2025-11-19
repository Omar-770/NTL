#pragma once
#include "optimiser.h"
#include <omp.h>
#include <iostream>
#include <chrono>

namespace NTL
{
	class NTL_opt;
	struct NTL_opt_setup;
	struct NTL_opt_result;


	struct NTL_opt_setup : public opt_setup
	{
		double Z0{ 0 };
		double er{ 0 };
		double d{ 0 };
		double Zs{ 0 };
		std::vector<double> Zl{};
		std::vector<double> freqs{};
		int K{ 0 };
		double m_Z_min{ 0 };
		double m_Z_max{ 0 };
	};

	struct NTL_opt_result : public opt_result
	{
		NTL ntl;
	};


	class NTL_opt : public opt
	{
	public:
		NTL_opt(const NTL_opt_setup& setup);

		NTL_opt_result optimise(console mode = console::inactive);

		NTL_opt_result optimise_d(double resolution, console mode = console::inactive);
		NTL_opt_result optimise_d(console mode = console::inactive);
		

	private:
		double m_Z0;
		double m_er;
		double m_d;
		double m_Zs;
		std::vector<double> m_Zl;
		std::vector<double> m_freqs;
		int m_K;
		double m_Z_min;
		double m_Z_max;

	private:
		//nlopt functions
		double min_objective(const std::vector<double>& Cn) const override;
		void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const override; 
	    void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const override;
		double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const override;
	};
}