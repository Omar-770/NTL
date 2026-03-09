#pragma once

#include <string>
#include <filesystem>
#include "models/wpd.h"
#include "models/ntl.h"

namespace NTL::latex
{
	bool ntl(const NTL& ntl, const std::string& folder, const std::string& file);
	bool wpd(const WPD::WPD& wpd, const NTL& out2, const NTL& out3, const std::string& folder, const std::string& file);
	bool sparams(const WPD::WPD& wpd, const NTL& out2, const NTL& out3,
		double Zref, double f_min, double f_max, double f_step, const std::vector<double>& freqs,
		const std::string& folder, const std::string& file);

	bool table_data(const WPD::WPD& wpd, const NTL& out2, const NTL& out3,
		const std::string& folder, const std::string& file);
}