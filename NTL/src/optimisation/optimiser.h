#pragma once

#include <vector>
#include <nlopt.hpp>
#include <omp.h>
#include "models/ntl.h"
#include "common/helpers.h"

namespace NTL
{
	class opt;
	struct opt_setup;
	struct opt_result;
	enum class console;

	struct opt_setup
	{
		int N{ 0 };
		std::vector<double> lb{};
		std::vector<double> ub{};
		std::vector<double> toll_bounds{};
		std::vector<double> toll_z{};
		double GBL_MAX{ 0 };
		double LCL_MAX{ 0 };
		double accepted_error{ 0 };
		int max_attempts{ 10'000 };
	};

	struct opt_result
	{
		std::vector<double> optimised_cn;
		double final_error;
	};

	enum class console
	{
		active = true, inactive = false
	};

	class opt
	{
	public:
		opt(const opt_setup& setup);
			

	protected:
		opt_result optimiser(console mode = console::inactive);
		double m_accepted_error;
		int m_max_attempts;

	private:
		int m_N;
		std::vector<double> m_lb;
		std::vector<double> m_ub;
		std::vector<double> m_toll_bounds;
		std::vector<double> m_toll_z;
		double m_GBL_MAX;
		double m_LCL_MAX;
		
		


	private:
		//nlopt functions
		virtual double min_objective(const std::vector<double>& Cn) const = 0;
		virtual void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const = 0;

	};
}