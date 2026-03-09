#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <complex>
#include <algorithm>
#include "qt_plot.h"
#include "models/ntl.h"
#include "simulation/NTL_sim.h"
#include "common/file_handler.h"
#include "common/enums.h"
#include "common/helpers.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);

		std::string file_name;
		std::cout << "File Name: ";
		std::cin >> file_name;
		std::cout << "\n\n\n";

		std::complex<double> Zs{};
		std::complex<double> Zl{};

		double f_min{};
		double f_max{};

		std::cout << "Source Impedance: ";
		std::cin >> Zs;
		std::cout << "Load Impedance: ";
		std::cin >> Zl;

		std::cout << "Minumim Frequency: ";
		std::cin >> f_min;

		std::cout << "Maximum Frequency: ";
		std::cin >> f_max;

		//build

		NTL::NTL ntl = fh::file_to_ntl(file_name);

		//simulate
		NTL::sim nsim(f_min, f_max);
		nsim.w_h_profile(ntl);
		nsim.z_profile(ntl);

		nsim.zin(ntl, Zl).magnitude();
		nsim.zin(ntl, Zl).phase();

		nsim.sparam(ntl, { Zs }, { Zl }).all();

		nsim.merge("NTL");

		app.exec();
	}
	catch (const std::exception& e)
	{
		std::cerr << "\nCRITICAL ERROR: " << e.what() << std::endl;
		std::cerr << "Press Enter to exit..." << std::endl;
		std::cin.get();
	}

	return 0;
}