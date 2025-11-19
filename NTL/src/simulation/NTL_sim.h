#pragma once

#include <utility>
#include "qt_plot.h"
#include "models/ntl.h"

namespace NTL
{

	class NTL_sim
	{
	public:
		NTL_sim()
		{

		}

		NTL_sim(const std::vector<NTL>& ntl) :m_ntl(ntl)
		{

		}

		void add(const NTL ntl) { m_ntl.push_back(ntl); } //to do: add a version add(ntl1, ntl2,...)

		void set_f_sweep(double f_min, double f_max, double f_step = 1e7)
		{
			m_fmin = f_min; m_fmax = f_max; m_fstep = f_step;
		}

		QMainWindow* merge(const char* title = "NTL_sim");

		QMainWindow* z_profile(int m, const char* title = "Impedance Z(z)", double step_size = 1e-4);
		QMainWindow* w_h_profile(int m, const char* title = "W/H(z)", double step_size = 1e-4);

		std::vector<QMainWindow*> s_matrix(int m, double Zs, double Zl, const char* title = "");
		std::vector<QMainWindow*> s_matrix(int m, double Zs, std::vector<double> Zl,
			std::vector<std::string> labels = {}, const char* title = "");

		QMainWindow* s_matrix(int m, int index, double Zs, double Zl, const char* title = "");
		QMainWindow* s_matrix(int m, int index, double Zs, std::vector<double> Zl,
			std::vector<std::string> labels = {}, const char* title = "");

		std::vector<QMainWindow*> get_windows() const { return m_windows; }

	private:
		std::vector<NTL> m_ntl;
		double m_fmin;
		double m_fmax;
		double m_fstep;

		qt_plot m_plotter;
		std::vector<QMainWindow*> m_windows;

	};
}