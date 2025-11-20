#pragma once

#include <utility>
#include "qt_plot.h"
#include "models/ntl.h"

namespace NTL
{
	class NTL_sim;


	class NTL_sim
	{
	public:
		NTL_sim() {};

		void set_f_sweep(double f_min, double f_max, double f_step = 1e7)
		{
			m_fmin = f_min; m_fmax = f_max; m_fstep = f_step;
		}

		QMainWindow* merge(const char* title = "NTL_sim");

		QMainWindow* z_profile(const NTL& ntl, const char* title = "Impedance Z(z)", double step_size = 1e-4);
		QMainWindow* w_h_profile(const NTL& ntl, const char* title = "W/H(z)", double step_size = 1e-4);

		std::vector<QMainWindow*> s_matrix(const NTL& ntl, double Zs, double Zl, const char* title = "");
		std::vector<QMainWindow*> s_matrix(const NTL& ntl, double Zs, std::vector<double> Zl,
			std::vector<std::string> labels = {}, const char* title = "");

		QMainWindow* s_matrix(const NTL& ntl, int index, double Zs, double Zl, const char* title = "");
		QMainWindow* s_matrix(const NTL& ntl
			, int index, double Zs, std::vector<double> Zl,
			std::vector<std::string> labels = {}, const char* title = "");

		std::vector<QMainWindow*> get_windows() const { return m_windows; }

	private:
		double m_fmin;
		double m_fmax;
		double m_fstep;

		qt_plot m_plotter;
		std::vector<QMainWindow*> m_windows;

	};
}