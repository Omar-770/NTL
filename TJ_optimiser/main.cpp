#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <omp.h>
#include <complex>
#include <algorithm>
#include <filesystem>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/tj.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/TJ_opt.h"
#include "simulation/NTL_sim.h"
#include "simulation/TJ_sim.h"
#include "common/file_handler.h"
#include "common/enums.h"
#include "common/helpers.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);

		//setup and run optimisation
		auto setup = fh::file_to_setup<TJ::opt_setup>("tj_setup");
		TJ::opt opt(setup);
		TJ::opt_result result = opt.optimise();

		//build the result
		TJ::TJ tj = result.tj;

		//simulate the result
		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		double f_min = 1e7;
		TJ::sim tsim(f_min, f_max * 1.2, 1e6, setup.freqs);

		tsim.w_h_profile(tj);

		tsim.z_profile(tj);

		std::array<std::complex<double>, 3> Zref;
		Zref[0] = { setup.Zref, 0 };
		Zref[1] = { setup.Zref, 0 };
		Zref[2] = { setup.Zref, 0 };
	
		tsim.sparams(tj, Zref, setup.K);
		tsim.merge("TJ");

		app.exec();

		//Saving

		std::cout << "\n\n\nSave Results? [y/n]: ";
		std::string save;
		std::cin >> save;
		if (save == "y" || save == "Y")
		{
			std::string folder;
			std::cout << "Folder Name: ";
			std::cin >> folder;

			std::cout << "Saving objects...\n";
			double H = 1.6e-3;
			//json
			std::filesystem::create_directory(folder);
			fh::tj_to_file(tj, folder + "/" + folder + "_tj");
			fh::ntl_to_file(tj.get_ntl2(), folder + "/" + folder + "_ntl2");
			fh::ntl_to_file(tj.get_ntl3(), folder + "/" + folder + "_ntl3");
			fh::setup_to_file<TJ::opt_setup>(setup, folder + "/" + folder + "_setup");

			std::cout << "\n\n\n=== Save Complete ===\n\n\n";
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << "\nCRITICAL ERROR: " << e.what() << std::endl;
		std::cerr << "Press Enter to exit..." << std::endl;
		std::cin.get();
	}

	return 0;
}