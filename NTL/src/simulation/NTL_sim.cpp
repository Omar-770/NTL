#include "NTL_sim.h"
#include <string>
namespace NTL
{
	QMainWindow* NTL_sim::merge(const char* title)
	{
		return m_plotter.combine(m_windows, title);
	}

	QMainWindow* NTL_sim::z_profile(const NTL& ntl,const char* title, double step_size)
	{
		auto z_vec = ntl.get_Z_vec(step_size);
		QMainWindow* window = m_plotter.plot(z_vec, "Z", title);
		m_windows.push_back(window);
		return window;
	}

	QMainWindow* zin_wrapper::magnitude(mag mode)
	{
		if (m_Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, Zin simulation " + std::string(m_title)));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin >= m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(m_title)));

		std::vector<std::pair<double, double>> Zin;
		std::vector<std::complex<double>> z = data();
		Zin.reserve(z.size());
		auto it = z.begin();

		for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
		{
			Zin.emplace_back(f, std::abs(*it));
			it++;
		}

		if (m_title[0] == '\0')
			m_title = mode == mag::abs ? "|Zin(f)|" : "|Zin(f)| dB";

		QMainWindow* window = m_sim.m_plotter.plot(Zin, "Zin", m_title);
		m_sim.m_windows.push_back(window);
		return window;
	}

	QMainWindow* zin_wrapper::phase(enum class phase mode)
	{
		if (m_Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, Zin simulation " + std::string(m_title)));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin >= m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(m_title)));

		std::vector<std::pair<double, double>> Zin;
		std::vector<std::complex<double>> z = data();
		Zin.reserve(z.size());
		auto it = z.begin();

		for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
		{
			Zin.emplace_back(f, (mode == phase::deg ? 180/M_PI : 1 ) * std::atan2(it->imag(), it->real()));
			it++;
		}

		QMainWindow* window = m_sim.m_plotter.plot(Zin, "Zin", m_title);
		m_sim.m_windows.push_back(window);
		return window;
	}

	std::vector<std::complex<double>> zin_wrapper::data()
	{
		std::vector<std::complex<double>> temp;
		temp.reserve((m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1);
		for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
		{
			temp.emplace_back(m_ntl.Zin(m_Zl, f));
		}

		return temp;
	}

	QMainWindow* NTL_sim::w_h_profile(const NTL& ntl,const char* title, double step_size)
	{
		auto w_h_vec = ntl.get_w_h_vec(step_size);
		QMainWindow* window = m_plotter.plot_mirror(w_h_vec, title);
		m_windows.push_back(window);
		return window;
	}

	QMainWindow* s_wrapper_2::magnitude(mag mode)
	{
		for(auto& z : m_Zl)
			if (z < 1e-6)
				throw(std::invalid_argument("Invalid terminal impedances, Zin simulation " + std::string(m_title)));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin >= m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(m_title)));

		std::vector<std::vector<std::pair<double, double>>> S;
		std::vector<std::complex<double>> temp;
		std::vector<const char*> S_label;
		std::vector<std::string> _temp;

		double points = (m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1;
		
		S.resize(m_Zl.size()); _temp.resize(m_Zl.size()); S_label.resize(m_Zl.size());

		for (int i = 0; i < m_Zl.size(); i++)
		{
			S[i].reserve(points);

			if (m_labels.empty())
				_temp[i] = "Z" + std::to_string(i + 1);
			else
				_temp[i] = m_labels[i];

			S_label[i] = _temp[i].c_str();			
		}

		temp.reserve(points);
		
		for (int i = 0; i < m_Zl.size(); i++)
		{
			temp = data(m_Zl[i]);
			auto it = temp.begin();
			for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
			{
				 if(mode == mag::dB)
					 S[i].emplace_back(f, 20 * std::log10(std::abs(*it)));
				 else
					 S[i].emplace_back(f, std::abs(*it));
				 it++;
			}
		}

		int first_index = m_index / 10;
		int second_index = m_index % 10;
		std::string title_temp;
		if (m_title[0] == '\0')
		{
			title_temp = "S" + std::to_string(first_index) + std::to_string(second_index);
			m_title = title_temp.c_str();
		}

		QMainWindow* window = m_sim.m_plotter.plot(S, S_label, m_title);
		m_sim.m_windows.push_back(window);
		return window;
	}

	QMainWindow* s_wrapper_2::phase(enum class phase mode)
	{
		for (auto& z : m_Zl)
			if (z < 1e-6)
				throw(std::invalid_argument("Invalid terminal impedances, Zin simulation " + std::string(m_title)));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin >= m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(m_title)));

		std::vector<std::vector<std::pair<double, double>>> S;
		std::vector<std::complex<double>> temp;
		std::vector<const char*> S_label;
		std::vector<std::string> _temp;

		double points = (m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1;

		S.resize(m_Zl.size()); _temp.resize(m_Zl.size()); S_label.resize(m_Zl.size());

		for (int i = 0; i < m_Zl.size(); i++)
		{
			S[i].reserve(points);

			if (m_labels.empty())
				_temp[i] = "Z" + std::to_string(i + 1);
			else
				_temp[i] = m_labels[i];

			S_label[i] = _temp[i].c_str();
		}

		temp.reserve(points);

		for (int i = 0; i < m_Zl.size(); i++)
		{
			temp = data(m_Zl[i]);
			auto it = temp.begin();
			for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
			{
				if (mode == phase::rad)
					S[i].emplace_back(f, std::atan2(it->imag(), it->real()));
				else
					S[i].emplace_back(f, 180 / M_PI * std::atan2(it->imag(), it->real()));
				it++;
			}
		}

		int first_index = m_index / 10;
		int second_index = m_index % 10;
		std::string title_temp;
		if (m_title[0] == '\0')
		{
			title_temp = "S" + std::to_string(first_index) + std::to_string(second_index);
			m_title = title_temp.c_str();
		}

		QMainWindow* window = m_sim.m_plotter.plot(S, S_label, m_title);
		m_sim.m_windows.push_back(window);
		return window;
	}

	std::vector<std::complex<double>> s_wrapper_2::data(double Zl)
	{	
		int first_index = m_index / 10;
		int second_index = m_index % 10;
	
		std::vector<std::complex<double>> temp;
		temp.reserve((m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1);
		for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
		{
			temp.emplace_back(m_ntl.S_matrix(f, m_Zs, Zl)(first_index - 1, second_index - 1));
		}

		return temp;
	}

	std::vector<QMainWindow*> s_wrapper_1::all(mag mode)
	{
		// 1. Guard Clauses (Safety First)
		for (auto& z : m_Zl)
			if (z < 1e-6)
				throw(std::invalid_argument("Invalid terminal impedances, S-Param simulation"));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin >= m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep in S-Param simulation"));


		std::vector<std::vector<std::pair<double, double>>> vS11, vS12, vS21, vS22;

		size_t num_loads = m_Zl.size();
		double points = (m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1;

		vS11.resize(num_loads); vS12.resize(num_loads);
		vS21.resize(num_loads); vS22.resize(num_loads);

		std::vector<std::string> _str_labels(num_loads);
		std::vector<const char*> char_labels(num_loads);

		for (size_t i = 0; i < num_loads; i++)
		{
			vS11[i].reserve(points); vS12[i].reserve(points);
			vS21[i].reserve(points); vS22[i].reserve(points);

			if (m_labels.empty())
				_str_labels[i] = "Z" + std::to_string(i + 1);
			else
				_str_labels[i] = m_labels[i];

			char_labels[i] = _str_labels[i].c_str();
		}

		for (size_t i = 0; i < num_loads; i++)
		{
			double Zl_current = m_Zl[i];

			for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
			{
				auto S = m_ntl.S_matrix(f, m_Zs, Zl_current);

				double v11 = std::abs(S(0, 0));
				double v12 = std::abs(S(0, 1));
				double v21 = std::abs(S(1, 0));
				double v22 = std::abs(S(1, 1));

				if (mode == mag::dB)
				{
					v11 = 20 * std::log10(v11);
					v12 = 20 * std::log10(v12);
					v21 = 20 * std::log10(v21);
					v22 = 20 * std::log10(v22);
				}

				vS11[i].emplace_back(f, v11);
				vS12[i].emplace_back(f, v12);
				vS21[i].emplace_back(f, v21);
				vS22[i].emplace_back(f, v22);
			}
		}

		std::vector<QMainWindow*> windows;
		QMainWindow* window;
		window = m_sim.m_plotter.plot(vS11, char_labels, "S11");
		windows.push_back(window);
		m_sim.m_windows.push_back(window);
		window = m_sim.m_plotter.plot(vS12, char_labels, "S12");
		windows.push_back(window);
		m_sim.m_windows.push_back(window);
		window = m_sim.m_plotter.plot(vS21, char_labels, "S21");
		m_sim.m_windows.push_back(window);
		windows.push_back(window);
		window = m_sim.m_plotter.plot(vS22, char_labels, "S22");
		m_sim.m_windows.push_back(window);
		windows.push_back(window);

		return windows;
	}
}