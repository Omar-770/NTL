#pragma once
#include "optimiser.h"
#include <omp.h>

namespace NTL
{
	struct NTL_opt_setup : public opt_setup
	{
		double Z0{ 0 };
		double er{ 0 };
		double d{ 0 };
		double Zs{ 0 };
		double Zl{ 0 };
		std::vector<double> freqs{};
		int K{ 0 };
		double m_Z_min{ 0 };
		double m_Z_max{ 0 };
	};


	class NTL_opt : public opt
	{
	public:
		NTL_opt(const NTL_opt_setup& setup) : opt(setup),
			m_Z0(setup.Z0), m_er(setup.er), m_d(setup.d), m_Zs(setup.Zs), m_Zl(setup.Zl),
			m_freqs(setup.freqs), m_K(setup.K), m_Z_min(setup.m_Z_min), m_Z_max(setup.m_Z_max)
		{

		}

		NTL optimise(bool output = false);
		NTL optimise(NTL& ntl, bool output = false);

	private:
		double m_Z0;
		double m_er;
		double m_d;
		double m_Zs;
		double m_Zl;
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