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

		//Setup Arms

		NTL::NTL ntl1 = fh::file_to_ntl("NTL1");
		NTL::NTL ntl2 = fh::file_to_ntl("NTL2");

		//Create WPD
		WPD::WPD wpd(ntl1, ntl2, 100);

		//Simulation
		WPD::sim wpdsim(1e7, 2.2e9, 1e6);


		wpdsim.w_h_profile(wpd);
		wpdsim.sparams(wpd, { 50, 50, 50});

		wpdsim.merge("WPD");

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