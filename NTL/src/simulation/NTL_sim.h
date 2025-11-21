#pragma once

#include <utility>
#include "qt_plot.h"
#include "models/ntl.h"

namespace NTL
{
	class NTL_sim;
	struct zin_wrapper;

	enum class phase
	{
		rad, deg
	};
	
	struct zin_wrapper
	{
		zin_wrapper(NTL_sim& sim, const NTL& ntl, double Zl)
			: m_sim(sim), m_ntl(ntl), m_Zl(Zl) {};
		QMainWindow* magnitude(const char* title = "Impedance |Zin(f)|");
		QMainWindow* phase(enum class phase mode = phase::deg, const char* title = "Impedance phase <Zin(f)");

	private:
		NTL_sim& m_sim;
		NTL m_ntl; double m_Zl;
		std::vector<std::complex<double>> data();
		
	};

	class NTL_sim
	{
	public:
		NTL_sim() : m_fmin(0), m_fmax(0), m_fstep(0) {};
		NTL_sim(double f_min, double f_max, double f_step = 1e7)
			: m_fmin(f_min), m_fmax(f_max), m_fstep(f_step)
		{

		}

		void set_f_sweep(double f_min, double f_max, double f_step = 1e7)
		{
			m_fmin = f_min; m_fmax = f_max; m_fstep = f_step;
		}

		QMainWindow* merge(const char* title = "NTL_sim");

		QMainWindow* z_profile(const NTL& ntl, const char* title = "Impedance Z(z)", double step_size = 1e-4);
		zin_wrapper zin(const NTL& ntl, double Zl);

		QMainWindow* w_h_profile(const NTL& ntl, const char* title = "W/H(z)", double step_size = 1e-4);

		std::vector<QMainWindow*> s_matrix(const NTL& ntl, double Zs, double Zl, const char* title = "");
		std::vector<QMainWindow*> s_matrix(const NTL& ntl, double Zs, std::vector<double> Zl,
			std::vector<std::string> labels = {}, const char* title = "");

		QMainWindow* s_matrix(const NTL& ntl, int index, double Zs, double Zl, const char* title = "");
		QMainWindow* s_matrix(const NTL& ntl
			, int index, double Zs, std::vector<double> Zl,
			std::vector<std::string> labels = {}, const char* title = "");

		std::vector<QMainWindow*> get_windows() const { return m_windows; }


		friend struct zin_wrapper;

	private:
		double m_fmin;
		double m_fmax;
		double m_fstep;

		qt_plot m_plotter;
		std::vector<QMainWindow*> m_windows;

	};
}