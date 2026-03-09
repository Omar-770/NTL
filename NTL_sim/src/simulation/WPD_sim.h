#pragma once

#include <utility>
#include "qt_plot.h"
#include "models/ntl.h"
#include "simulation/NTL_sim.h"
#include "models/wpd.h"
#include "common/enums.h"

namespace WPD
{
	class sim;

	class sim
	{
	public:
		sim() : m_fmin(0), m_fmax(0), m_fstep(0) {};

		sim(double f_min, double f_max, double f_step = 1e6, const std::vector<double>& freqs_target = {})
			: m_fmin(f_min), m_fmax(f_max), m_fstep(f_step)
		{
			set_target_f(freqs_target);
		}

		void set_f_sweep(double f_min, double f_max, double f_step = 1e6)
		{
			m_fmin = f_min; m_fmax = f_max; m_fstep = f_step;
		}

		void set_target_f(const std::vector<double>& f)
		{
			m_freqs = f;
			std::sort(m_freqs.begin(), m_freqs.end());
		}

		void add_window(QMainWindow* window)
		{
			m_windows.push_back(window);
		}

		QMainWindow* merge(const char* title = "WPD_sim");

		std::vector<QMainWindow*> z_profile(const WPD& wpd);
		std::vector<QMainWindow*> w_h_profile(const WPD& wpd);

		std::vector<QMainWindow*> sparams(const WPD& wpd, std::array<std::complex<double>, 3> Zl, int K = 50);
		std::vector<QMainWindow*> sparams(const WPD& wpd, double Zref, NTL::NTL& output2, NTL::NTL& output3, int K = 50);

		std::vector<QMainWindow*> get_windows() const { return m_windows; }

	private:
		double m_fmin;
		double m_fmax;
		double m_fstep;

		qt_plot m_plotter;
		std::vector<QMainWindow*> m_windows;
		NTL::sim NTL_sim;

		std::vector<double> m_freqs;
	};
}