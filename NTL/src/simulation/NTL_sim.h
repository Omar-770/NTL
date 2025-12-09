#pragma once

#include <utility>
#include "qt_plot.h"
#include "models/ntl.h"
#include "common/enums.h"

namespace NTL
{
	class sim;
	struct zin_wrapper;
	struct s_wrapper_1;
	struct s_wrapper_2;
	
	struct zin_wrapper
	{
		zin_wrapper(sim& sim, const NTL& ntl, double Zl, const char* title)
			: m_sim(sim), m_ntl(ntl), m_Zl(Zl), m_title(title) {};
		QMainWindow* magnitude(mag mode = mag::abs);
		QMainWindow* phase(enum class phase mode = phase::deg);

	private:
		sim& m_sim;
		const NTL& m_ntl; double m_Zl;
		const char* m_title;
		std::vector<std::complex<double>> data();
		
	};

	struct s_wrapper_2
	{
		s_wrapper_2(sim& sim, const NTL& ntl, double Zs, const std::vector<double>& Zl, 
			const std::vector<std::string>& labels, const char* title, int index)
			: m_sim(sim), m_ntl(ntl), m_Zs(Zs), m_Zl(Zl), m_labels(labels), m_title(title), m_index(index) {
		};

		QMainWindow* magnitude(mag mode = mag::dB);
		QMainWindow* phase(enum class phase mode = phase::deg);
	private:
		sim& m_sim;
		const NTL& m_ntl; double m_Zs; const std::vector<double>& m_Zl;
		const std::vector<std::string>& m_labels; const char* m_title;
		int m_index;
		std::vector<std::complex<double>> data(double Zl);
	};

	struct s_wrapper_1
	{
		s_wrapper_1(sim& sim, const NTL& ntl, double Zs, const std::vector<double>& Zl,
			const std::vector<std::string>& labels, const char* title)
			: m_sim(sim), m_ntl(ntl), m_Zs(Zs), m_Zl(Zl), m_labels(labels), m_title(title) {};

		s_wrapper_2 S11() { return s_wrapper_2(m_sim, m_ntl, m_Zs, m_Zl, m_labels, m_title, 11); }
		s_wrapper_2 S12() { return s_wrapper_2(m_sim, m_ntl, m_Zs, m_Zl, m_labels, m_title, 12); }
		s_wrapper_2 S21() { return s_wrapper_2(m_sim, m_ntl, m_Zs, m_Zl, m_labels, m_title, 21); }
		s_wrapper_2 S22() { return s_wrapper_2(m_sim, m_ntl, m_Zs, m_Zl, m_labels, m_title, 22); }
		std::vector<QMainWindow*> all(mag mode = mag::dB);
	private:
		sim& m_sim;
		const NTL& m_ntl; double m_Zs; const std::vector<double>& m_Zl;
		const std::vector<std::string>& m_labels; const char* m_title;
	};


	class sim
	{
	public:
		sim() : m_fmin(0), m_fmax(0), m_fstep(0) {};
		sim(double f_min, double f_max, double f_step = 1e6, const std::vector<double>& freqs_target = {})
			: m_fmin(f_min), m_fmax(f_max), m_fstep(f_step)
		{
			set_target_f(freqs_target);
		}

		void set_f_sweep(double f_min, double f_max, double f_step = 1e7)
		{
			m_fmin = f_min; m_fmax = f_max; m_fstep = f_step;
		}

		void set_target_f(const std::vector<double>& f)
		{
			m_freqs = f;
			std::sort(m_freqs.begin(), m_freqs.end());
		}

		QMainWindow* merge(const char* title = "NTL_sim");

		QMainWindow* z_profile(const NTL& ntl, const char* title = "Impedance Z(z)", double step_size = 1e-4);
		zin_wrapper zin(const NTL& ntl, double Zl, const char* title = "") { return zin_wrapper(*this, ntl, Zl, title); }

		QMainWindow* w_h_profile(const NTL& ntl, const char* title = "W/H(z)", double step_size = 1e-4);

		s_wrapper_1 sparam(const NTL& ntl, double Zs, const std::vector<double>& Zl,
			std::vector<std::string> labels = {}, const char* title = "") { return s_wrapper_1(*this, ntl, Zs, Zl, labels, title); }

		std::vector<QMainWindow*> get_windows() const { return m_windows; }

		
		friend struct zin_wrapper;
		friend struct s_wrapper_1;
		friend struct s_wrapper_2;

	private:
		double m_fmin;
		double m_fmax;
		double m_fstep;

		qt_plot m_plotter;
		std::vector<QMainWindow*> m_windows;
		std::vector<double> m_freqs;

	};
}