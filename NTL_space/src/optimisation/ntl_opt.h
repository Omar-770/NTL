#pragma once

#include <vector>
#include <nlopt.hpp>
#include <omp.h>
#include "models/ntl.h"
#include "common/helpers.h"

namespace NTL
{
	class NTL_opt;

	struct NTL_opt_setup
	{
		double Z0;
		double er;
		double d;
		double Zs;
		double Zl;
		std::vector<double> freqs;
		int K;

		int N;
		double Z_min;
		double Z_max;

		double toll_bounds;
		double toll_z;
		double GBL_MAX;
		double LCL_MAX;
		double accepted_error;
	};

	class NTL_opt
	{
	public:
		NTL_opt(const NTL_opt_setup& setup)
			: m_Z0(setup.Z0), m_er(setup.er), m_d(setup.d),
			m_Zs(setup.Zs), m_Zl(setup.Zl), m_freqs(setup.freqs),
			m_K(setup.K), m_N(setup.N), m_Z_min(setup.Z_min),
			m_Z_max(setup.Z_max), m_GBL_MAX(setup.GBL_MAX), m_LCL_MAX(setup.LCL_MAX),
			m_accepted_error(setup.accepted_error)
		{
			m_lb = std::vector<double>(m_N, -1);
			m_ub = std::vector<double>(m_N, 1);
			m_toll_bounds = std::vector<double>(2, setup.toll_bounds);
			m_toll_z = std::vector<double>(m_K, setup.toll_z);
		}

		NTL optimise(bool output = false);


	private:
		double m_Z0;
		double m_er;
		double m_d;
		double m_Zs;
		double m_Zl;
		std::vector<double> m_freqs;
		int m_K;
		int m_N;
		double m_Z_min;
		double m_Z_max;

		std::vector<double> m_lb;
		std::vector<double> m_ub;
		std::vector<double> m_toll_bounds;
		std::vector<double> m_toll_z;
		double m_GBL_MAX;
		double m_LCL_MAX;
		double m_accepted_error;

	public:
		NTL_opt_setup get_setup() const
		{
			return NTL_opt_setup{
				m_Z0,m_er,m_d,m_Zs,m_Zl,m_freqs,m_K,m_N,m_Z_min,
				m_Z_max,m_toll_bounds[0],m_toll_z[0],m_GBL_MAX,
				m_LCL_MAX,m_accepted_error
			};
		}

	private:
		//nlopt functions

		double min_objective(const std::vector<double>& Cn) const;
		void equality_constraints(unsigned m, double* res, unsigned n, const double* Cn) const;
		void inequality_constraints_Zmax(unsigned m, double* res, unsigned n, const double* Cn) const;
		void inequality_constraints_Zmin(unsigned m, double* res, unsigned n, const double* Cn) const;

		double objective_with_fd_gradient(const std::vector<double>& Cn, std::vector<double>& grad, void* data) const;

	};
}