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

	zin_wrapper NTL_sim::zin(const NTL& ntl, double Zl)
	{
		return zin_wrapper(*this, ntl, Zl);
	}

	QMainWindow* zin_wrapper::magnitude(const char* title)
	{
		if (m_Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin > m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(title)));

		std::vector<std::pair<double, double>> Zin;

		double points = (m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1;
		Zin.reserve(points);

		for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
		{
			Zin.emplace_back(f, std::abs(m_ntl.Zin(m_Zl, f)));
		}

		QMainWindow* window = m_sim.m_plotter.plot(Zin, "Zin", title);
		m_sim.m_windows.push_back(window);
		return window;
	}

	QMainWindow* zin_wrapper::phase(enum class phase mode,const char* title)
	{
		if (m_Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));

		if (!m_sim.m_fmin || !m_sim.m_fmax || !m_sim.m_fstep || m_sim.m_fmin > m_sim.m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(title)));

		std::vector<std::pair<double, double>> Zin;

		double points = (m_sim.m_fmax - m_sim.m_fmin) / m_sim.m_fstep + 1;
		Zin.reserve(points);

		for (double f = m_sim.m_fmin; f < m_sim.m_fmax; f += m_sim.m_fstep)
		{
			std::complex<double> temp = m_ntl.Zin(m_Zl, f);
			Zin.emplace_back(f, ((mode == phase::deg) ? 180 / M_PI : 1) * std::atan2(temp.imag(), temp.real()));
		}

		QMainWindow* window = m_sim.m_plotter.plot(Zin, "Zin", title);
		m_sim.m_windows.push_back(window);
		return window;
	}

	QMainWindow* NTL_sim::w_h_profile(const NTL& ntl,const char* title, double step_size)
	{
		auto w_h_vec = ntl.get_w_h_vec(step_size);
		QMainWindow* window = m_plotter.plot_mirror(w_h_vec, title);
		m_windows.push_back(window);
		return window;
	}

	std::vector<QMainWindow*> NTL_sim::s_matrix(const NTL& ntl, double Zs, double Zl, const char* title)
	{
		if (Zs < 1e-6 || Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));

		if (!m_fmin || !m_fmax || !m_fstep || m_fmin > m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(title)));

		std::vector<std::pair<double, double>> S11, S12, S21, S22;

		double points = (m_fmax - m_fmin) / m_fstep + 1;
		S11.reserve(points);
		S12.reserve(points);
		S21.reserve(points);
		S22.reserve(points);

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			matrix2x2cd S_matrix = ntl.S_matrix(f, Zs, Zl);

			S11.emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 0))));
			S12.emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 1))));
			S21.emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 0))));
			S22.emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 1))));
		}

		std::vector<QMainWindow*> window_vec;

		window_vec.push_back(m_plotter.plot(S11, "S11", ("S11 " + std::string(title)).c_str()));
		window_vec.push_back(m_plotter.plot(S12, "S12", ("S12 " + std::string(title)).c_str()));
		window_vec.push_back(m_plotter.plot(S21, "S21", ("S21 " + std::string(title)).c_str()));
		window_vec.push_back(m_plotter.plot(S22, "S22", ("S22 " + std::string(title)).c_str()));

		for (int i = 0; i < window_vec.size(); i++)
			m_windows.push_back(window_vec[i]);


		return window_vec;
	}

	std::vector<QMainWindow*> NTL_sim::s_matrix(const NTL& ntl, double Zs, std::vector<double> Zl, std::vector<std::string> labels, const char* title)
	{
		if (Zs < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));
		for (auto& Z : Zl)
			if (Z < 1e-6)
				throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));
		if (!labels.empty() && labels.size() != Zl.size())
			throw(std::invalid_argument("Invalid number of labels, S_matrix simulation " + std::string(title)));
		if (!m_fmin || !m_fmax || !m_fstep || m_fmin > m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(title)));

		std::vector<std::vector<std::pair<double, double>>> S11, S12, S21, S22;
		std::vector<const char*> S11_label, S12_label, S21_label, S22_label;
		std::vector<std::string> temp;
		S11.resize(Zl.size()); S11_label.resize(Zl.size()); temp.resize(Zl.size());
		S12.resize(Zl.size()); S12_label.resize(Zl.size());
		S21.resize(Zl.size()); S21_label.resize(Zl.size());
		S22.resize(Zl.size()); S22_label.resize(Zl.size());

		double points = (m_fmax - m_fmin) / m_fstep + 1;

		for (int i = 0; i < Zl.size(); i++)
		{
			S11[i].reserve(points);
			S12[i].reserve(points);
			S21[i].reserve(points);
			S22[i].reserve(points);

			if (labels.empty())
				temp[i] = "Z" + std::to_string(i + 1);
			else
				temp[i] = labels[i];
			S11_label[i] = temp[i].c_str();
			S12_label[i] = temp[i].c_str();
			S21_label[i] = temp[i].c_str();
			S22_label[i] = temp[i].c_str();
		}

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			for (int i = 0; i < Zl.size(); i++)
			{
				matrix2x2cd S_matrix = ntl.S_matrix(f, Zs, Zl[i]);

				S11[i].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 0))));
				S12[i].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 1))));
				S21[i].emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 0))));
				S22[i].emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 1))));
			}
		}

		std::vector<QMainWindow*> window_vec;

		window_vec.push_back(m_plotter.plot(S11, S11_label, ("S11 " + std::string(title)).c_str()));
		window_vec.push_back(m_plotter.plot(S12, S12_label, ("S12 " + std::string(title)).c_str()));
		window_vec.push_back(m_plotter.plot(S21, S21_label, ("S21 " + std::string(title)).c_str()));
		window_vec.push_back(m_plotter.plot(S22, S22_label, ("S22 " + std::string(title)).c_str()));

		for (int i = 0; i < window_vec.size(); i++)
			m_windows.push_back(window_vec[i]);

		return window_vec;
	}

	QMainWindow* NTL_sim::s_matrix(const NTL& ntl, int index, double Zs, double Zl, const char* title)
	{
		if (Zs < 1e-6 || Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation" + std::string(title)));

		int first_index = index / 10;
		int second_index = index % 10;

		if (first_index < 1 || first_index > 2 ||
			second_index < 1 || second_index > 2)
			throw(std::invalid_argument("Invalid indices, S_matrix simulation " + std::string(title)));
		if (!m_fmin || !m_fmax || !m_fstep || m_fmin > m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(title)));

		std::vector<std::pair<double, double>> S;

		double points = (m_fmax - m_fmin) / m_fstep + 1;
		S.reserve(points);

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			matrix2x2cd S_matrix = ntl.S_matrix(f, Zs, Zl);

			S.emplace_back(f, 20 * std::log10(std::abs(S_matrix(first_index - 1, second_index - 1))));
		}

		std::string title_temp;
		if (title[0] == '\0')
		{
			title_temp = "S" + std::to_string(first_index) + std::to_string(second_index);
			title = title_temp.c_str();
		}
		QMainWindow* window = m_plotter.plot(S, ("S" + std::to_string(first_index) + std::to_string(second_index)).c_str(), title);
		m_windows.push_back(window);
		return window;
	}

	QMainWindow* NTL_sim::s_matrix(const NTL& ntl, int index, double Zs, std::vector<double> Zl, std::vector<std::string> labels, const char* title)
	{
		if (Zs < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));
		for (auto& Z : Zl)
			if (Z < 1e-6)
				throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation " + std::string(title)));
		if (!labels.empty() && labels.size() != Zl.size())
			throw(std::invalid_argument("Invalid number of labels, S_matrix simulation " + std::string(title)));
		if (!m_fmin || !m_fmax || !m_fstep || m_fmin > m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep " + std::string(title)));

		int first_index = index / 10;
		int second_index = index % 10;

		if (first_index < 1 || first_index > 2 ||
			second_index < 1 || second_index > 2)
			throw(std::invalid_argument("Invalid indices, S_matrix simulation " + std::string(title)));

		std::vector<std::vector<std::pair<double, double>>> S;
		std::vector<const char*> S_label;
		std::vector<std::string> temp;
		S.resize(Zl.size()); S_label.resize(Zl.size()); temp.resize(Zl.size());

		double points = (m_fmax - m_fmin) / m_fstep + 1;

		for (int i = 0; i < Zl.size(); i++)
		{
			S[i].reserve(points);

			if (labels.empty())
				temp[i] = "Z" + std::to_string(i + 1);
			else
				temp[i] = labels[i];
			S_label[i] = temp[i].c_str();

		}

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			for (int i = 0; i < Zl.size(); i++)
			{
				matrix2x2cd S_matrix = ntl.S_matrix(f, Zs, Zl[i]);

				S[i].emplace_back(f, 20 * std::log10(std::abs(S_matrix(first_index - 1, second_index - 1))));
			}
		}

		std::string title_temp;
		if (title[0] == '\0')
		{
			title_temp = "S" + std::to_string(first_index) + std::to_string(second_index);
			title = title_temp.c_str();
		}
		QMainWindow* window = m_plotter.plot(S, S_label, title);

		m_windows.push_back(window);
		return window;
	}
}