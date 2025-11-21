#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "qt_plot.h"
#include "models/ntl.h"
#include "models/wpd.h"
#include "optimisation/NTL_opt.h"
#include "optimisation/WPD_opt.h"
#include "simulation/NTL_sim.h"
#include "simulation/WPD_sim.h"
#include "common/file_handler.h"

namespace fh = NTL::fh;

int main(int argc, char* argv[])
{
	try
	{
		QApplication app(argc, argv);
		auto start_time = std::chrono::high_resolution_clock::now();
		auto setup1 = fh::file_to_setup<NTL::NTL_opt_setup>("setup1");
		auto setup2 = fh::file_to_setup<NTL::NTL_opt_setup>("setup2");

		NTL::NTL ntl1 = fh::file_to_ntl("NTL1");
		NTL::NTL ntl2 = fh::file_to_ntl("NTL2");

		NTL::NTL_sim sim(1e6, 2.2e9, 1e6);

		sim.w_h_profile(ntl1, "NTL1");
		sim.w_h_profile(ntl2, "NTL2");

		sim.zin(ntl1, 50).magnitude();
		sim.zin(ntl2, 50).magnitude();

		sim.sparam(ntl1, setup1.Zs, setup1.Zl).all();
		//sim.sparam(ntl1, setup1.Zs, setup1.Zl).S11().magnitude();
		//sim.sparam(ntl1, setup1.Zs, setup1.Zl).S11().phase();
	
		//sim.sparam(ntl2, setup2.Zs, setup2.Zl).S11().magnitude();
		//sim.sparam(ntl2, setup2.Zs, setup2.Zl).S11().phase();

		sim.merge("NTL");
		auto end_time = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = end_time - start_time;
		std::cout << "\n\n\n*** Execution finished in " << elapsed.count() / 60.0 << " minutes" << std::endl;
		app.exec();
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