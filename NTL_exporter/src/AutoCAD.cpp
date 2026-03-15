#include <fstream>
#include <iostream>
#include "AutoCAD.h"

void scr(const NTL::NTL& ntl, double substrate_height, const std::string& filename, double step)
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