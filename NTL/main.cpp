#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <omp.h>
#include <complex>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/wpd.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/WPD_opt.h"
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

		WPD::opt opt(setup);
		WPD::opt_result result = opt.optimise();

		WPD::WPD wpd = result.wpd;
		NTL::NTL out2 = result.output2;
		NTL::NTL out3 = result.output3;

		WPD::sim sim(1e7, 2.2e9, 1e6, setup.freqs);
		NTL::sim out_sim(1e7, 2.2e9, 1e6);
		sim.w_h_profile(wpd);
		sim.add_window(out_sim.w_h_profile(out2, "Out2"));
		sim.add_window(out_sim.w_h_profile(out3, "Out3"));
		sim.sparams(wpd, setup.Zref, out2, out3);	

		sim.merge("WPD");

		fh::wpd_to_file(wpd, "single_band_equal_wpd/wpd");
		fh::ntl_to_file(out2, "single_band_equal_wpd/out2");
		fh::ntl_to_file(out3, "single_band_equal_wpd/out3");
		fh::setup_to_file<WPD::opt_setup>(setup, "single_band_equal_wpd/setup");
		

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