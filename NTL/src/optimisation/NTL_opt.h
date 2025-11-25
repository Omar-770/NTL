#pragma once

#include "optimiser.h"
#include <omp.h>
#include <iostream>
#include <nlohmann/json.hpp>
#include <chrono>

namespace NTL
{
	class NTL_opt;
	struct NTL_opt_setup;
	struct NTL_opt_result;
	

	struct NTL_opt_setup : public opt_setup
	{
		NTL_opt_setup() {};
		NTL_opt_setup(const NTL_opt_setup& setup) : opt_setup(setup) {};
		NTL_opt_setup(const nlohmann::json& j);
		double Z0{ 0 };
		double er{ 0 };
		double d{ 0 };
		double Zs{ 0 };
		std::vector<double> Zl{};
		std::vector<double> freqs{};
		int K{ 0 };
		double Z_min{ 0 };
		double Z_max{ 0 };

		nlohmann::json get_json() const override;
	};

	struct NTL_opt_result : public opt_result
	{
		NTL ntl;
	};


	class NTL_opt : public opt
	{
	public:
		NTL_opt(const NTL_opt_setup& setup);

		NTL_opt_result optimise(console mode = console::active);

		NTL_opt_result optimise_d(double resolution, console mode = console::active);
		NTL_opt_result optimise_d(console mode = console::active);
		

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