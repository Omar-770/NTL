#include "models/ntl.h"
#include "models/wpd.h"
#include "models/tj.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/WPD_opt.h"
#include "optimisation/TJ_opt.h"

inline double verify_result(const NTL::opt_setup setup, const NTL::opt_result& result)
{
	double matching_error{};
	auto Zl = setup.Zl;
	auto Zs = setup.Zs;

	if (setup.Zl.size() == 1)
	{
		Zl.resize(setup.freqs.size());
		for (int i = 1; i < Zl.size(); i++)
			Zl[i] = Zl[0];
	}
	if (setup.Zs.size() == 1)
	{
		Zs.resize(setup.freqs.size());
		for (int i = 1; i < Zs.size(); i++)
			Zs[i] = Zs[0];
	}

	//matching error
	for (int i = 0; i < setup.freqs.size(); i++)
	{
		auto S = result.ntl.S_matrix(setup.freqs[i], Zs[i], Zl[i], setup.K);
		double matching = 20 * std::log10(std::max(std::abs(S(0, 0)), 1e-15));
		if (matching >= -10)
		{
			matching_error += std::pow(matching + 10, 2);
		}
	}

	matching_error = std::sqrt(matching_error);

	std::cout << "Matching Error: " << matching_error << (matching_error == 0 ? " -- Perfect\n" : "\n");

	return matching_error;
}

inline double verify_result(const WPD::opt_setup setup, const WPD::opt_result& result)
{
	double total_error{};
	double matching_error{};
	double split_error{};
	double isolation_error{};

	auto split = setup.split;
	if (setup.split.size() == 1)
	{
		split.resize(setup.freqs.size());
		for (int i = 1; i < split.size(); i++)
			split[i] = split[0];
	}

	for (int i = 0; i < setup.freqs.size(); i++)
	{
		auto S = result.wpd.system_S_matrix(setup.Zref, result.output2, result.output3, setup.freqs[i], setup.K);

		//matching error
		for(int j = 0; j < 3; j++)
		{
			double matching = 20 * std::log10(std::max(std::abs(S(j, j)), 1e-15));
			if (matching >= -10)
			{
				matching_error += std::pow(matching + 10, 2);
			}
		}

		//split error
		double s21_dB_target = 10 * std::log10(1.0 / (split[i] + 1));
		double s31_dB_target = 10 * std::log10(split[i] / (split[i] + 1));
		
		double s21_dB = 20 * std::log10(std::max(std::abs(S(1, 0)), 1e-15));
		double s31_dB = 20 * std::log10(std::max(std::abs(S(2, 0)), 1e-15));

		split_error += std::pow(s21_dB - s21_dB_target, 2);
		split_error += std::pow(s31_dB - s31_dB_target, 2);

		//isolation error
		double isolation = 20 * std::log10(std::max(std::abs(S(2, 1)), 1e-15));
		if (isolation >= -20)
		{
			isolation_error += std::pow(isolation + 10, 2);
		}
	}

	total_error = matching_error + split_error + isolation_error;
	matching_error = std::sqrt(matching_error);
	split_error = std::sqrt(split_error);
	isolation_error = std::sqrt(isolation_error);
	total_error = std::sqrt(total_error);	

	std::cout << "Matching Error: " << matching_error << (matching_error == 0 ? " -- Perfect\n" : "\n");
	std::cout << "Split Error: " << split_error << (split_error == 0 ? " -- Perfect\n" : "\n");
	std::cout << "Isolation Error: " << isolation_error << (isolation_error == 0 ? " -- Perfect\n" : "\n");
	std::cout << "Total Error: " << total_error << (total_error == 0 ? " -- Perfect!!!\n" : "\n");
	return total_error;
}

inline double verify_result(const TJ::opt_setup setup, const TJ::opt_result& result)
{
	double total_error{};
	double matching_error{};
	double split_error{};

	std::array<std::complex<double>, 3> Zref;
	Zref[0] = { setup.Zref, 0 };
	Zref[1] = { setup.Zref, 0 };
	Zref[2] = { setup.Zref, 0 };


	auto split = setup.split;
	if (setup.split.size() == 1)
	{
		split.resize(setup.freqs.size());
		for (int i = 1; i < split.size(); i++)
			split[i] = split[0];
	}

	for (int i = 0; i < setup.freqs.size(); i++)
	{
		auto S = result.tj.S_matrix(setup.freqs[i], Zref, setup.K);

		//matching error
		for (int j = 0; j < 3; j++)
		{
			double matching = 20 * std::log10(std::max(std::abs(S(j, j)), 1e-15));
			if (matching >= -10)
			{
				matching_error += std::pow(matching + 10, 2);
			}
		}

		//split error
		double s21_dB_target = 10 * std::log10(1.0 / (split[i] + 1));
		double s31_dB_target = 10 * std::log10(split[i] / (split[i] + 1));

		double s21_dB = 20 * std::log10(std::max(std::abs(S(1, 0)), 1e-15));
		double s31_dB = 20 * std::log10(std::max(std::abs(S(2, 0)), 1e-15));

		split_error += std::pow(s21_dB - s21_dB_target, 2);
		split_error += std::pow(s31_dB - s31_dB_target, 2);

	}

	total_error = matching_error + split_error;
	matching_error = std::sqrt(matching_error);
	split_error = std::sqrt(split_error);
	total_error = std::sqrt(total_error);

	std::cout << "Matching Error: " << matching_error << (matching_error == 0 ? " -- Perfect\n" : "\n");
	std::cout << "Split Error: " << split_error << (split_error == 0 ? " -- Perfect\n" : "\n");
	std::cout << "Total Error: " << total_error << (total_error == 0 ? " -- Perfect!!!\n" : "\n");
	return total_error;
}