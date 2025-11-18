#include "NTL_sim.h"
#include <string>
namespace NTL
{
	void NTL_sim::merge(const char* title)
	{
		m_plotter.combine(m_windows, title);
	}

	void NTL_sim::z_profile(const char* title, double step_size)
	{
		auto z_vec = m_ntl.get_Z_vec(step_size);
		m_windows.push_back(m_plotter.plot(z_vec, title));
	}

	void NTL_sim::w_h_profile(const char* title, double step_size)
	{
		auto w_h_vec = m_ntl.get_w_h_vec(step_size);

		m_windows.push_back(m_plotter.plot_mirror(w_h_vec, title));
	}

	void NTL_sim::s_matrix(double Zs, double Zl, const char* title)
	{
		if (Zs < 0 || Zl < 0 || 0 < Zs < 1e-6 || 0 < Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminsl impedances, S_matrix simulation" + std::string(title)));

		std::vector<std::pair<double, double>> S11, S12, S21, S22;
		
		double points = (m_fmax - m_fmin) / m_fstep + 1;
		S11.reserve(points);
		S12.reserve(points);
		S21.reserve(points);
		S22.reserve(points);

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			matrix2x2cd S_matrix = m_ntl.S_matrix(f, Zs, Zl);

			S11.emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 0))));
			S12.emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 1))));
			S21.emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 0))));
			S22.emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 1))));
		}

		m_windows.push_back(m_plotter.plot(S11, ("S11" + std::string(title)).c_str()));
		m_windows.push_back(m_plotter.plot(S12, ("S12" + std::string(title)).c_str()));
		m_windows.push_back(m_plotter.plot(S21, ("S21" + std::string(title)).c_str()));
		m_windows.push_back(m_plotter.plot(S22, ("S22" + std::string(title)).c_str()));
	}

	void NTL_sim::s_matrix(int index, double Zs, double Zl, const char* title)
	{
		if (Zs < 0 || Zl < 0 || 0 < Zs < 1e-6 || 0 < Zl < 1e-6)
			throw(std::invalid_argument("Invalid terminal impedances, S_matrix simulation" + std::string(title)));

		int first_index = index / 10;
		int second_index = index % 10;

		if (first_index < 1 || first_index > 2 ||
			second_index < 1 || second_index > 2)
			throw(std::invalid_argument("Invalid indices, S_matrix simulation " + std::string(title)));

		std::vector<std::pair<double, double>> S;
		
		double points = (m_fmax - m_fmin) / m_fstep + 1;
		S.reserve(points);

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			matrix2x2cd S_matrix = m_ntl.S_matrix(f, Zs, Zl);

			S.emplace_back(f, 20 * std::log10(std::abs(S_matrix(first_index - 1, second_index - 1))));
		}

		m_windows.push_back(m_plotter.plot(S, ("S" + std::to_string(index) + std::string(title)).c_str()));

	}
}