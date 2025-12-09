#pragma once

#include "optimiser.h"
#include <omp.h>
#include <iostream>
#include <nlohmann/json.hpp>
#include <chrono>
#include <algorithm>

namespace NTL
{
	class opt;
	struct opt_setup;
	struct opt_result;
	

	struct opt_setup : public optimiser_setup
	{
		opt_setup() {};
		opt_setup(const opt_setup& setup);
		opt_setup(const optimiser_setup& setup);
		opt_setup(const nlohmann::json& j);
		double Z0{ 0 };
		double er{ 0 };
		double d{ 0 };
		int M{ 0 };
		double Zs{ 0 };
		std::vector<double> Zl{};
		std::vector<double> freqs{};
		int K{ 0 };
		double Z_min{ 0 };
		double Z_max{ 0 };
		double Z_at_0{ 0 };
		double Z_at_d{ 0 };
	
		nlohmann::json get_json() const override;
	};

	struct opt_result : public optimiser_result
	{
		NTL ntl;
	};


	class opt : public optimiser
	{
	public:
		opt(const opt_setup& setup);

		opt_result optimise(console mode = console::active);

		opt_result optimise_d(double resolution, console mode = console::active);
		opt_result optimise_d(console mode = console::active);
		

	private:
		double m_Z0;
		double m_er;
		double m_d;
		int m_M;
		double m_Zs;
		std::vector<double> m_Zl;
		std::vector<double> m_freqs;
		int m_K;
		double m_Z_min;
		double m_Z_max;
		double m_Z_at_0;
		double m_Z_at_d;

	private:
		//nlopt functions
		double min_objective(const std::vector<double>& Cn) const override;
		void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const override;
		void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const override; 
	    void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const override;
		double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const override;
	};
}