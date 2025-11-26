#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <omp.h>
#include <algorithm>
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
		auto start_time = std::chrono::high_resolution_clock::now();

		
		NTL::NTL ntl = fh::file_to_ntl("NTL1");
		WPD::WPD wpd(ntl, ntl, 200);
		WPD::sim sim(1e7, 2.2e9, 1e6);
		NTL::sim simn(1e7, 2.2e9, 1e6);

		std::array<std::complex<double>, 3> Zl;

		simn.sparam(ntl, 100, { 50 }, {}, "100 to 50").S11().magnitude();
		simn.sparam(ntl, 50, { 100 }, {}, "50 to 100").S11().magnitude();

		sim.w_h_profile(wpd);
		sim.sparams(wpd, {50, 100, 100});
		

		sim.merge("WPD");
		simn.merge("NTL");

		auto end_time = std::chrono::high_resolution_clock::now();
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