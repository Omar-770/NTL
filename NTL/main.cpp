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
		setup1.d = 80e-3;
		setup1.Zl = { 100 };

		NTL::NTL_opt opt(setup1);		
		NTL::NTL ntl = opt.optimise().ntl;
		NTL::NTL_sim sim(1e6, 2.2e9, 1e6);

		sim.w_h_profile(ntl);
		sim.zin(ntl, 50).magnitude();
		sim.zin(ntl, 50).phase();
		sim.sparam(ntl, setup1.Zs, setup1.Zl).S11().magnitude();
		
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