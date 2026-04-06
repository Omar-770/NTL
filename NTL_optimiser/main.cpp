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
#include "optimisation/NTL_opt.h"
#include "simulation/NTL_sim.h"
#include "common/file_handler.h"
#include "common/enums.h"
#include "common/helpers.h"
#include "validation.h"
#include "verification.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);

		//setup and run optimisation
		auto setup = fh::file_to_setup<NTL::opt_setup>("ntl_setup");
		NTL::opt opt(setup);
		NTL::opt_result result = opt.optimise();		

		//build the result
		NTL::NTL ntl = result.ntl;

		//validate the result
		std::cout << "\nValidating Result...\n";
		if (validate_result(setup, result))
			std::cout << "Valid Result...\n";
		else
			std::cout << "Invalid Result...\n";

		std::cout << "\n===============================\n";

		//verify the result
		std::cout << "\nVerifying Result...\n";
		verify_result(setup, result);

		std::cout << "\n===============================\n";

		//simulate the result
		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		double f_min = 1e7;
		NTL::sim nsim(f_min, f_max * 1.2, 1e6, setup.freqs);

		nsim.w_h_profile(ntl);
		nsim.z_profile(ntl);

		nsim.sparam(ntl, setup.Zs, setup.Zl).all();

		nsim.merge("NTL");

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
			fh::ntl_to_file(ntl, folder + "/" + folder + "_ntl");
			fh::setup_to_file<NTL::opt_setup>(setup, folder + "/" + folder + "_setup");

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