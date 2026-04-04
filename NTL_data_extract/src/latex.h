#pragma once

#include <string>
#include <filesystem>
#include "models/wpd.h"
#include "models/ntl.h"
#include "models/tj.h"

bool ntl_profiles(const NTL::NTL& ntl, const std::string& folder, const std::string& file);
bool ntl_table_data(const NTL::NTL& ntl, const std::string& folder, const std::string& file);
bool wpd_profiles(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3, const std::string& folder, const std::string& file);
bool wpd_sparams(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3,
			double Zref, double f_min, double f_max, double f_step, const std::vector<double>& freqs,
		    const std::string& folder, const std::string& file);

bool wpd_table_data(const WPD::WPD& wpd, const NTL::NTL& out2, const NTL::NTL& out3,
	const std::string& folder, const std::string& file);

bool tj_profiles(const TJ::TJ& tj, const std::string& folder, const std::string& file);
bool tj_sparams(const TJ::TJ& tj, double Zref, double f_min, double f_max, double f_step, const std::vector<double>& freqs,
		    const std::string& folder, const std::string& file);
bool tj_table_data(const TJ::TJ& tj, const std::string& folder, const std::string& file);