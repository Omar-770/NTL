#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <omp.h>
#include <complex>
#include <algorithm>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/wpd.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/WPD_opt.h"
#include "optimisation/WPD_opt_2.h"
#include "simulation/NTL_sim.h"
#include "simulation/WPD_sim.h"
#include "common/file_handler.h"
#include "common/enums.h"
#include "common/helpers.h"
#include "common/latex.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);
		
		//setup and run optimisation
		auto setup = fh::file_to_setup<WPD::opt_setup>("wpd_setup1");
		WPD::opt_2 opt(setup);
		WPD::opt_result result = opt.optimise();

		//build the result
		NTL::NTL out2 = result.output2;
		NTL::NTL out3 = result.output3;
		WPD::WPD wpd = result.wpd;

		//simulate the result
		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		double f_min = 1e7;
		WPD::sim wsim(f_min, f_max * 1.2, 1e6, setup.freqs);
		NTL::sim nsim;

		wsim.w_h_profile(wpd);
		wsim.add_window(nsim.w_h_profile(out2, "Out2"));
		wsim.add_window(nsim.w_h_profile(out3, "Out3"));

		wsim.z_profile(wpd);
		wsim.add_window(nsim.z_profile(out2, "Out2-Z"));
		wsim.add_window(nsim.z_profile(out3, "Out3-Z"));

		wsim.sparams(wpd, setup.Zref, out2, out3);
		wsim.merge("WPD");

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
			//Latex

			std::cout << "Saving to CSVs...\n";

			if (NTL::latex::wpd(wpd, out2, out3, folder, folder) &&
				NTL::latex::sparams(wpd, out2, out3, setup.Zref, f_min, f_max, 1e6, setup.freqs, folder, folder) &&
				NTL::latex::table_data(wpd, out2, out3, folder, folder))
				std::cout << "CSVs save success\n";
			else
				std::cerr << "CSVs save failure\n";

			std::cout << "Saving objects...\n";
			double H = 1.6e-3;
			//json
			fh::wpd_to_file(wpd, folder + "/" + folder + "_wpd");
			fh::ntl_to_file(out2, folder + "/" + folder + "_out2");
			fh::ntl_to_file(out3, folder + "/" + folder + "_out3");
			fh::setup_to_file<WPD::opt_setup>(setup, folder + "/" + folder + "_setup");
			//AutoCAD Script
			fh::export_geometry_scr(out2, H, folder + "/" + folder + "_out2");
			fh::export_geometry_scr(out3, H, folder + "/" + folder + "_out3");
			fh::export_geometry_scr(wpd.get_ntl2(), H, folder + "/" + folder + "_ntl2");
			fh::export_geometry_scr(wpd.get_ntl3(), H, folder + "/" + folder + "_ntl3");

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