#include "file_handler.h"
#include "helpers.h"


namespace NTL::fh
{
	std::fstream json_to_file(const json& j, const std::string& name)
	{
		std::fstream stream;
		stream.open(name + ".json", std::ios_base::out);

		if (stream)
		{
			stream << j.dump();
		}
		else
			throw(std::runtime_error("Couldn't open " + name + " for writing"));

		stream.close();

		return stream;
	}

	json file_to_json(const std::string& name)
	{
		std::fstream stream;
		stream.open(name + ".json", std::ios_base::in);
		json j;

		if (stream)
		{
			j = json::parse(stream);
		}
		else
			throw(std::runtime_error("Couldn't open " + name + " for reading"));

		return j;
	}

	json ntl_to_json(const NTL& ntl, const std::string& readme)
	{
		if (ntl.get_Cn().empty())
			throw(std::logic_error("Attempted to write an empty NTL into a json object"));

		return {
			{"json_type", "NTL"},
			{"Z0", ntl.get_Z0()},
			{"e_r", ntl.get_er()},
			{"d", ntl.get_d()},
			{"Cn", ntl.get_Cn()},
			{"M", ntl.get_M()},
			{"readme", readme}
		};
	}

	NTL json_to_ntl(const json& j)
	{
		if (j.at("json_type") != "NTL")
			throw(std::logic_error("Attempted to read an NTL from a different json object"));

		return NTL(j.at("Z0"), j.at("e_r"), j.at("d"), j.at("Cn"), j.at("M"));
	}

	json wpd_to_json(const WPD::WPD& wpd, const std::string& readme)
	{
		if (wpd.get_ntl2().get_Cn().empty() || wpd.get_ntl3().get_Cn().empty() || !wpd.get_R())
			throw(std::logic_error("Attempted to write an empty NTL into a json object"));

		json j1 = {
			{"json_type", "WPD"},
			{"Z0", wpd.get_Z0()},
			{"e_r", wpd.get_er()},
			{"R", wpd.get_R()},
			{"readme", readme}
		};

		json j2 = ntl_to_json(wpd.get_ntl2());
		j2.erase("json_type");
		j2.erase("readme");
		j2.erase("Z0");
		j2.erase("e_r");
		rename_json_element<double>(j2, "d", "d-2");
		rename_json_element<std::vector<double>>(j2, "Cn", "Cn2");
		rename_json_element<double>(j2, "M", "M-2");

		json j3 = ntl_to_json(wpd.get_ntl3());
		j3.erase("json_type");
		j3.erase("readme");
		j3.erase("Z0");
		j3.erase("e_r");
		rename_json_element<double>(j3, "d", "d-3");
		rename_json_element<std::vector<double>>(j3, "Cn", "Cn3");
		rename_json_element<double>(j3, "M", "M-3");

		json j_total = j1;
		j_total.merge_patch(j2);
		j_total.merge_patch(j3);

		return j_total;
	}

	WPD::WPD json_to_wpd(const json& j)
	{
		if (j.at("json_type") != "WPD")
			throw(std::logic_error("Attempted to read a WPD from a different json object"));

		return WPD::WPD(NTL(j.at("Z0"), j.at("e_r"), j.at("d-2"), j.at("Cn2"), j.at("M-2")),
			NTL(j.at("Z0"), j.at("e_r"), j.at("d-3"), j.at("Cn3"), j.at("M-3")),
			j.at("R"));
	}

	json tj_to_json(const TJ::TJ& tj, const std::string& readme)
	{
		if (tj.get_ntl2().get_Cn().empty() || tj.get_ntl3().get_Cn().empty())
			throw(std::logic_error("Attempted to write an empty NTL into a json object"));

		json j1 = {
			{"json_type", "TJ"},
			{"Z0", tj.get_Z0()},
			{"e_r", tj.get_er()},
			{"readme", readme}
		};

		json j2 = ntl_to_json(tj.get_ntl2());
		j2.erase("json_type");
		j2.erase("readme");
		j2.erase("Z0");
		j2.erase("e_r");
		rename_json_element<double>(j2, "d", "d-2");
		rename_json_element<std::vector<double>>(j2, "Cn", "Cn2");
		rename_json_element<double>(j2, "M", "M-2");

		json j3 = ntl_to_json(tj.get_ntl3());
		j3.erase("json_type");
		j3.erase("readme");
		j3.erase("Z0");
		j3.erase("e_r");
		rename_json_element<double>(j3, "d", "d-3");
		rename_json_element<std::vector<double>>(j3, "Cn", "Cn3");
		rename_json_element<double>(j3, "M", "M-3");

		json j_total = j1;
		j_total.merge_patch(j2);
		j_total.merge_patch(j3);

		return j_total;;
	}

	TJ::TJ json_to_tj(const json& j)
	{
		if (j.at("json_type") != "TJ")
			throw(std::logic_error("Attempted to read a WPD from a different json object"));

		return TJ::TJ(NTL(j.at("Z0"), j.at("e_r"), j.at("d-2"), j.at("Cn2"), j.at("M-2")),
			NTL(j.at("Z0"), j.at("e_r"), j.at("d-3"), j.at("Cn3"), j.at("M-3")));
	}

	void export_geometry_csv(const NTL& ntl, double substrate_height, const std::string& filename, double step)
	{
		// Get the (z, W/H) points from the NTL model
		auto profile = ntl.get_w_h_vec(step);

		std::fstream stream;
		stream.open(filename + ".csv", std::ios_base::out);

		if (!stream.is_open())
			throw std::runtime_error("Could not open file for CSV export");

		// Header for clarity
		stream << "X_mm,Width_mm\n";

		for (const auto& point : profile)
		{
			double z_meters = point.first;
			double w_h_ratio = point.second;

			// Convert to Millimeters for HFSS convenience
			double x_mm = z_meters * 1000.0;
			double w_mm = w_h_ratio * substrate_height * 1000.0;

			stream << x_mm << "," << w_mm << "\n";
		}

		stream.close();
		std::cout << "Geometry exported to " << filename << ".csv" << std::endl;
	}

	void export_geometry_scr(const NTL& ntl, double substrate_height, const std::string& filename, double step)
	{
		// 1. Get the profile data (Position z, Width/Height ratio)
		auto profile = ntl.get_w_h_vec(step);

		std::fstream stream;
		stream.open(filename + ".scr", std::ios_base::out);

		if (!stream.is_open())
			throw std::runtime_error("Could not open file for AutoCAD export");

		// 2. Setup AutoCAD Environment (Optional but recommended)
		stream << "_OSMODE 0\n"; // Turn off Object Snap to avoid drawing errors
		stream << "_UCS _World\n"; // Ensure we are in World Coordinates

		// 3. Start the Polyline command
		stream << "_PLINE\n";

		// 4. Draw Top Edge (Forward)
		// We convert meters to millimeters for CAD (standard convention)
		for (const auto& point : profile)
		{
			double x_mm = point.first * 1000.0;
			double w_mm = point.second * substrate_height * 1000.0;

			// AutoCAD Format: X,Y
			stream << x_mm << "," << (w_mm / 2.0) << "\n";
		}

		// 5. Draw Bottom Edge (Backward) to close the loop
		// Iterate in reverse
		for (auto it = profile.rbegin(); it != profile.rend(); ++it)
		{
			double x_mm = it->first * 1000.0;
			double w_mm = it->second * substrate_height * 1000.0;

			stream << x_mm << "," << (-w_mm / 2.0) << "\n";
		}

		// 6. Close the Polyline and Zoom Extents
		stream << "_CLOSE\n"; // 'c' closes the polyline back to the start
		stream << "_ZOOM _E\n"; // Zoom to show the whole part

		stream.close();
		std::cout << "AutoCAD script exported to " << filename << ".scr" << std::endl;
	}

	void bend_geometry_scr(const NTL& ntl, double substrate_height, const std::string& filename, double angle, double step)
	{
		if (angle <= 0 || angle > 360)
			throw std::invalid_argument("Angle must be between 0 and 360");

		std::vector<std::pair<double, double>> profile = ntl.get_w_h_vec(step);
		size_t size = profile.size();

		substrate_height *= 1000.0;

		std::fstream stream;
		stream.open(filename + "_" + std::to_string(int(angle)) + "_bend" + ".scr", std::ios_base::out);

		if (!stream.is_open())
			throw std::runtime_error("Could not open file for AutoCAD export");

		stream << "_OSMODE 0\n";
		stream << "_UCS _World\n";
		stream << "_PLINE\n";

		double Theta = (angle / 180.0 * M_PI);
		double R = ntl.get_d() * 1000.0 / Theta;

		double step_theta = Theta / (size - 1);

		for (size_t i = 0; i < size; i++)
		{
			double theta = double(i) * step_theta;
			if (theta > Theta)
			{
				theta = Theta;
			}

			double x_mm = (R - profile[i].second * substrate_height / 2.0) * std::cos(theta);
			double y_mm = (R - profile[i].second * substrate_height / 2.0) * std::sin(theta);

			stream << x_mm << "," << y_mm << "\n";
		}

		for (int i = size - 1; i >= 0; i--)
		{
			double theta = double(i) * step_theta;
			if (theta < 0)
			{
				theta = 0;
			}

			double x_mm = (R + profile[i].second * substrate_height / 2.0) * std::cos(theta);
			double y_mm = (R + profile[i].second * substrate_height / 2.0) * std::sin(theta);

			stream << x_mm << "," << y_mm << "\n";
		}

		stream << "_CLOSE\n";
		stream << "_REGION\n_L\n\n";
		stream << "_ZOOM _E\n";

		stream.close();
		std::cout << "AutoCAD script exported to " << (filename + "_" + std::to_string(int(angle)) + "_bend") << ".scr" << std::endl;
		std::cout << std::endl;
	}

	std::fstream ntl_to_file(const NTL& ntl, const std::string& name, const std::string& readme)
	{
		return json_to_file(ntl_to_json(ntl, readme), name);
	}

	NTL file_to_ntl(const std::string& name)
	{
		return json_to_ntl(file_to_json(name));
	}

	std::fstream wpd_to_file(const WPD::WPD& wpd, const std::string& name, const std::string& readme)
	{
		return json_to_file(wpd_to_json(wpd, readme), name);
	}
	WPD::WPD file_to_wpd(const std::string& name)
	{
		return json_to_wpd(file_to_json(name));
	}
	std::fstream tj_to_file(const TJ::TJ& tj, const std::string& name, const std::string& readme)
	{
		return json_to_file(tj_to_json(tj, readme), name);
	}
	TJ::TJ file_to_tj(const std::string& name)
	{
		return json_to_tj(file_to_json(name));
	}
}