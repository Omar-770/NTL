#pragma once


#include <nlohmann/json.hpp>
#include <fstream>
#include "models/ntl.h"
#include "models/wpd.h"
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

	//WPD to json
	json wpd_to_json(const WPD::WPD& ntl, const std::string& readme = std::string());
	WPD::WPD json_to_wpd(const json& j);

	//setup to json
	template<class T>
	json setup_to_json(const T& setup, const std::string& readme = std::string())
	{
		return setup.get_json();
	}

	template<class T>
	T json_to_setup(const json& j)
	{
		return T(j);
	}

	/// Shortcut functions
	//NTL to file
	std::fstream ntl_to_file(const NTL& ntl, const std::string& name, const std::string& readme = std::string());
	NTL file_to_ntl(const std::string& name);

	//WPD to file
	std::fstream wpd_to_file(const WPD::WPD& wpd, const std::string& name, const std::string& readme = std::string());
	WPD::WPD file_to_wpd(const std::string& name);

	//setup to file
	template<typename T>
	std::fstream setup_to_file(const T& setup, const std::string& name, const std::string& readme = std::string())
	{
		return json_to_file(setup_to_json<T>(setup, readme), name);
	}

	template<typename T>
	T file_to_setup(const std::string& name)
	{
		return json_to_setup<T>(file_to_json(name));
	}

	//Helpers
	template<typename T>
	void rename_json_element(json& j, const std::string& old_name, const std::string& new_name)
	{
		T val = j.at(old_name).get<T>();
		j.erase(old_name);
		j[new_name] = val;
	}

}