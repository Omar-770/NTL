#include "WPD_sim.h"

namespace WPD
{
	QMainWindow* sim::merge(const char* title)
	{
		return m_plotter.combine(m_windows, title);
	}

	std::vector<QMainWindow*> sim::z_profile(const WPD& wpd)
	{
		QMainWindow* window1 = NTL_sim.z_profile(wpd.get_ntl2(), "NTL2");
		QMainWindow* window2 = NTL_sim.z_profile(wpd.get_ntl3(), "NTL3");
		m_windows.push_back(window1);
		m_windows.push_back(window2);

		return { window1, window2 };
	}

	std::vector<QMainWindow*> sim::w_h_profile(const WPD& wpd)
	{
		QMainWindow* window1 = NTL_sim.w_h_profile(wpd.get_ntl2(), "NTL2");
		QMainWindow* window2 = NTL_sim.w_h_profile(wpd.get_ntl3(), "NTL3");
		m_windows.push_back(window1);
		m_windows.push_back(window2);

		return { window1, window2 };
	}

	std::vector<QMainWindow*> sim::sparams(const WPD& wpd, std::array<std::complex<double>, 3> Zl, int K)
	{
		for (auto& z : Zl)
			if (std::abs(z) < 1e-12)
				throw(std::invalid_argument("Invalid terminal impedances, Zin simulation "));
		if (!m_fmin || !m_fmax || !m_fstep || m_fmin >= m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep "));

		std::vector<std::vector<std::pair<double, double>>> S;
		S.resize(6);
		double points = (m_fmax - m_fmin) / m_fstep + 1;
		for (auto& vec : S)
			vec.reserve(points);

		auto log_freq = m_freqs.cbegin();
		auto log_freq_end = m_freqs.cend();

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			matrix3x3cd S_matrix = wpd.S_matrix(f, Zl, K);


			S[0].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 0)))); //S11
			S[1].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 1)))); //S12
			S[2].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 2)))); //S13
			S[3].emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 1)))); //S22
			S[4].emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 2)))); //S23
			S[5].emplace_back(f, 20 * std::log10(std::abs(S_matrix(2, 2)))); //S33	

			if (log_freq != log_freq_end && std::abs(f - *log_freq) < 1e-6)
			{
				std::cout << "[Frequency: " << f << "]\n";
				std::cout << "Matching:\n";
				std::cout << "\tS11: " << 20 * std::log10(std::abs(S_matrix(0, 0))) << '\n';
				std::cout << "\tS22: " << 20 * std::log10(std::abs(S_matrix(1, 1))) << '\n';
				std::cout << "\tS33: " << 20 * std::log10(std::abs(S_matrix(2, 2))) << '\n';
				std::cout << "Split:\n";
				std::cout << "\tS12: " << 20 * std::log10(std::abs(S_matrix(0, 1))) << '\n';
				std::cout << "\tS13: " << 20 * std::log10(std::abs(S_matrix(0, 2))) << '\n';
				std::cout << "Isolation:\n";
				std::cout << "\tS23: " << 20 * std::log10(std::abs(S_matrix(1, 2))) << "\n\n";
				log_freq++;
			}
		}

		std::vector<QMainWindow*> windows;
		QMainWindow* window;

		window = m_plotter.plot({ S[0], S[3], S[5] }, { "S11", "S22", "S33" }, "Matching");
		windows.push_back(window);
		m_windows.push_back(window);
		window = m_plotter.plot({ S[1], S[2] }, { "S12", "S13" }, "Split");
		windows.push_back(window);
		m_windows.push_back(window);
		window = m_plotter.plot(S[4], "S23", "Isolation");
		windows.push_back(window);
		m_windows.push_back(window);

		return windows;
	}
	std::vector<QMainWindow*> sim::sparams(const WPD& wpd, double Zref, NTL::NTL& output2, NTL::NTL& output3, int K)
	{
		if (Zref < 1e-12)
			throw(std::invalid_argument("Invalid terminal impedances, Zin simulation "));
		if (!m_fmin || !m_fmax || !m_fstep || m_fmin >= m_fmax)
			throw(std::invalid_argument("Invalid frequency sweep "));

		std::vector<std::vector<std::pair<double, double>>> S;
		S.resize(6);
		double points = (m_fmax - m_fmin) / m_fstep + 1;
		for (auto& vec : S)
			vec.reserve(points);

		auto log_freq = m_freqs.cbegin();
		auto log_freq_end = m_freqs.cend();

		for (double f = m_fmin; f < m_fmax; f += m_fstep)
		{
			matrix3x3cd S_matrix = wpd.system_S_matrix(Zref, output2, output3, f, K);


			S[0].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 0)))); //S11
			S[1].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 1)))); //S12
			S[2].emplace_back(f, 20 * std::log10(std::abs(S_matrix(0, 2)))); //S13
			S[3].emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 1)))); //S22
			S[4].emplace_back(f, 20 * std::log10(std::abs(S_matrix(1, 2)))); //S23
			S[5].emplace_back(f, 20 * std::log10(std::abs(S_matrix(2, 2)))); //S33							

			if (log_freq != log_freq_end && std::abs(f - *log_freq) < 1e-6)
			{
				std::cout << "[Frequency: " << f << "]\n";
				std::cout << "Matching:\n";
				std::cout << "\tS11: " << 20 * std::log10(std::abs(S_matrix(0, 0))) << '\n';
				std::cout << "\tS22: " << 20 * std::log10(std::abs(S_matrix(1, 1))) << '\n';
				std::cout << "\tS33: " << 20 * std::log10(std::abs(S_matrix(2, 2))) << '\n';
				std::cout << "Split:\n";
				std::cout << "\tS12: " << 20 * std::log10(std::abs(S_matrix(0, 1))) << '\n';
				std::cout << "\tS13: " << 20 * std::log10(std::abs(S_matrix(0, 2))) << '\n';
				std::cout << "Isolation:\n";
				std::cout << "\tS23: " << 20 * std::log10(std::abs(S_matrix(1, 2))) << "\n\n";
				log_freq++;
			}

		}

		std::vector<QMainWindow*> windows;
		QMainWindow* window;

		window = m_plotter.plot({ S[0], S[3], S[5] }, { "S11", "S22", "S33" }, "Matching");
		windows.push_back(window);
		m_windows.push_back(window);
		window = m_plotter.plot({ S[1], S[2] }, { "S12", "S13" }, "Split");
		windows.push_back(window);
		m_windows.push_back(window);
		window = m_plotter.plot(S[4], "S23", "Isolation");
		windows.push_back(window);
		m_windows.push_back(window);

		return windows;
	}
}