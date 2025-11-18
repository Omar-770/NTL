#pragma once
#include "optimiser.h"
#include <omp.h>

namespace NTL
{
	class NTL_opt;
	struct NTL_opt_setup;


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


	class NTL_opt : public opt
	{
	public:
		NTL_opt(const NTL_opt_setup& setup) : opt(setup),
			m_Z0(setup.Z0), m_er(setup.er), m_d(setup.d), m_Zs(setup.Zs), m_Zl(setup.Zl),
			m_freqs(setup.freqs), m_K(setup.K), m_Z_min(setup.m_Z_min), m_Z_max(setup.m_Z_max)
		{
			if (m_Zl.size() < m_freqs.size() && m_Zl.size() != 1)
				throw(std::invalid_argument("Number of load impedances & frequency points mismatch"));

			if (m_Zl.size() == 1)
			{
				m_Zl.resize(m_freqs.size());
				for (int i = 1; i < m_Zl.size(); i++)
					m_Zl[i] = m_Zl[0];
			}
		}

		NTL optimise(console mode = console::inactive);
		NTL optimise(NTL& ntl, console mode = console::inactive);

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