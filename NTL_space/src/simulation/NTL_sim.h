#pragma once

#include <utility>
#include "qt_plot.h"
#include "models/ntl.h"

namespace NTL
{
	class NTL_sim
	{
	public:
		NTL_sim(const NTL& ntl) : m_ntl(ntl)
		{

		}

		void set_f_sweep(double f_min, double f_max, double f_step = 1e7)
		{
			m_fmin = f_min; m_fmax = f_max; m_fstep = f_step;
		}

		void merge(const char* title = "NTL_sim");

		void z_profile(double step_size = 1e-4);
		void w_h_profile(double step_size = 1e-4);
		void s_matrix(double Zs, double Zl);


	private:
		NTL m_ntl;
		double m_fmin;
		double m_fmax;
		double m_fstep;

		qt_plot m_plotter;
		std::vector<QMainWindow*> m_windows;
	};
}