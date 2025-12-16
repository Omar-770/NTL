#include <iostream>
#include <Eigen/Dense>
#include <vector>
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

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);

		auto setup = fh::file_to_setup<WPD::opt_setup>("wpd_setup1");

		WPD::opt_2 opt(setup);
		WPD::opt_result result = opt.optimise();

		WPD::WPD wpd = result.wpd;
		NTL::NTL out2 = result.output2;
		NTL::NTL out3 = result.output3;

		double f_max = *std::max_element(setup.freqs.cbegin(), setup.freqs.cend());
		f_max *= 1.5;

		WPD::sim sim(1e7, f_max, 1e6, setup.freqs);
		NTL::sim out_sim(1e7, f_max, 1e6);

		sim.w_h_profile(wpd);
		sim.add_window(out_sim.w_h_profile(out2, "Out2"));
		sim.add_window(out_sim.w_h_profile(out3, "Out3"));
		sim.sparams(wpd, setup.Zref, out2, out3);	

		sim.merge("WPD");

		opt.print_result_logs();
		

		app.exec();

		std::cout << "\n\n\n";
		std::string save;
		do
		{
			std::cout << "Save results? [Y/N]: ";			
			std::cin >> save;
		} while (save != "y" && save != "Y" && save != "N" && save != "n");
		bool saved = false;
		std::string folder;
		if(save == "Y" || save == "y")
		{					
			do
			{
				try
				{					
					std::cout << "Enter folder name: ";
					std::cin >> folder;
					if (folder == "0") break;
					fh::wpd_to_file(wpd, folder + "/wpd");
					fh::ntl_to_file(out2, folder + "/out2");
					fh::ntl_to_file(out3, folder + "/out3");
					fh::setup_to_file<WPD::opt_setup>(setup, folder + "/setup");
					saved = true;
				}
				catch (std::exception& e)
				{
					std::cerr << "Saving ERROR: " << e.what() << ", check if file exists...\n";
				}
			} while (saved != true);			
		}

		if (saved)
			std::cout << "\n\nResults saved under " << folder << std::endl;
		else
			std::cout << "\n\nExiting pogram without saving...";
	}
	catch (const std::exception& e)
	{
		std::cerr << "\nCRITICAL ERROR: " << e.what() << std::endl;
		std::cerr << "Press Enter to exit..." << std::endl;
		std::cin.get();
	}

	return 0;
}