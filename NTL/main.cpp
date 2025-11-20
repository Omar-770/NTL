#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "qt_plot.h"
#include "models/ntl.h"
#include "optimisation/optimiser.h"
#include "optimisation/NTL_opt.h"
#include "simulation/NTL_sim.h"
#include "common/file_handler.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		auto start_time = std::chrono::high_resolution_clock::now();
		QApplication app(argc, argv);

		NTL::NTL_opt_setup setup1, setup2;												

		setup1 = dynamic_cast<NTL::NTL_opt_setup&>(*fh::file_to_setup("setup1"));
		setup2 = dynamic_cast<NTL::NTL_opt_setup&>(*fh::file_to_setup("setup2"));

		NTL::NTL ntl1, ntl2;
		NTL::NTL_opt opt1(setup1), opt2(setup2);

		std::cout << "Optimising first NTL\n\n";
		ntl1 = opt1.optimise_d(0.5e-3, NTL::console::active).ntl;
		std::cout << "\n\nOptimising second NTL\n\n";
		ntl2 = opt2.optimise_d(0.5e-3, NTL::console::active).ntl;

		NTL::NTL_sim sim;

		sim.w_h_profile(ntl1, "W/H(z) (1)");
		sim.w_h_profile(ntl2, "W/H(z) (2)");

		sim.set_f_sweep(1e7, 2.2e9);
		sim.s_matrix(ntl1, 11, setup1.Zs, setup1.Zl, { "50", "100", "150" }, "S11 (1)");
		sim.s_matrix(ntl2, 11, setup2.Zs, setup2.Zl, { "50", "25", "16.67" }, "S11 (2)");

		sim.merge("NTL");

		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;

		std::cout << "\n\n\n*** Execution finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		app.exec();

		fh::ntl_to_file(ntl1, "NTL1", "matches 50 to 50/100/150 at 0.5/1.3/2");
		fh::ntl_to_file(ntl2, "NTL2", "matches 50 to 50/25/16.67 at 0.5/1.3/2");
	}
	catch (std::invalid_argument& e)
	{
		std::cout << "Invalid argument exception: ";
		std::cout << e.what() << std::endl;
	}
	catch (std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}

	return 0;
}