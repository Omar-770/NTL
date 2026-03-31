#pragma once

#include <string>
#include <filesystem>
#include "models/wpd.h"
#include "models/ntl.h"

bool ntl_profiles(const NTL::NTL& ntl, const std::string& folder, const std::string& file);
bool wpd_profiles(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3, const std::string& folder, const std::string& file);
bool wpd_sparams(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3,
			double Zref, double f_min, double f_max, double f_step, const std::vector<double>& freqs,
		    const std::string& folder, const std::string& file);

bool wpd_table_data(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3,
	const std::string& folder, const std::string& file);
