#pragma once


#include <nlohmann/json.hpp>
#include <fstream>
#include "models/ntl.h"
#include "simulation/NTL_sim.h"
#include "optimisation/NTL_opt.h"


namespace NTL::fh
{
	using json = nlohmann::json;

	//json to file
	std::fstream json_to_file(const json& j, const std::string& name);
	json file_to_json(const std::string& name);

	//NTL to json
	json ntl_to_json(const NTL& ntl, const std::string& readme = std::string());
	NTL json_to_ntl(const json& j);

	//setup to json
	json setup_to_json(const opt_setup& setup, const std::string& readme = std::string());
	std::unique_ptr<opt_setup> json_to_setup(const json& j);

	/// Shortcut functions
	//NTL to file
	std::fstream ntl_to_file(const NTL& ntl, const std::string& name, const std::string& readme = std::string());
	NTL file_to_ntl(const std::string& name);
	//setup to file
	std::fstream setup_to_file(const opt_setup& setup, const std::string& name, const std::string& readme = std::string());
	std::unique_ptr<opt_setup> file_to_setup(const std::string& name);


}