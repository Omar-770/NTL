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
	};

	class opt
	{
	public:
		opt(const opt_setup& setup)
			: m_N(setup.N), m_lb(setup.lb), m_ub(setup.ub), m_toll_bounds(setup.toll_bounds), m_toll_z(setup.toll_z),
			m_GBL_MAX(setup.GBL_MAX), m_LCL_MAX(setup.LCL_MAX), m_accepted_error(setup.accepted_error)
		{

		}

	protected:
		std::vector<double> optimiser(bool output = false);

	private:
		int m_N;
		std::vector<double> m_lb;
		std::vector<double> m_ub;
		std::vector<double> m_toll_bounds;
		std::vector<double> m_toll_z;
		double m_GBL_MAX;
		double m_LCL_MAX;
		double m_accepted_error;


	private:
		//nlopt functions
		virtual double min_objective(const std::vector<double>& Cn) const = 0;
		virtual void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const = 0;
		virtual double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const = 0;

	};
}