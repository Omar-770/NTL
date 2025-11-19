#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include "qt_plot.h"
#include "models/ntl.h"
#include "optimisation/optimiser.h"
#include "optimisation/NTL_opt.h"
#include "simulation/NTL_sim.h"


int main(int argc, char* argv[])
{
	try
	{

		QApplication app(argc, argv);

		NTL::NTL_opt_setup setup1, setup2;
		setup1.Z0 = 50; setup1.Zs = 50; setup1.Zl = { 50, 100, 150 }; setup1.er = 4.6; setup1.d = 66e-3; setup1.freqs = { 0.5e9, 1.3e9, 2e9 };
		setup1.N = 10; setup1.K = 50; setup1.lb = std::vector<double>(setup1.N, -1);  setup1.ub = std::vector<double>(setup1.N, 1);
		setup1.toll_bounds = std::vector<double>(2, 1e-6); setup1.toll_z = std::vector<double>(setup1.K, 1e-6);
		setup1.GBL_MAX = 350e3; setup1.LCL_MAX = 50e3; setup1.accepted_error = 6e-3;
		setup1.m_Z_min = 0.1 * setup1.Z0; setup1.m_Z_max = 2.6 * setup1.Z0; setup1.max_attempts = 20;

		setup2 = setup1;
		setup2.Zl = { 50, 50.0 / 2, 50.0 / 3 };
		setup2.d = 64e-3;

		NTL::NTL ntl1, ntl2;
		NTL::NTL_opt opt1(setup1), opt2(setup2);
		

		std::cout << "Optimising first NTL\n\n";
		ntl1 = opt1.optimise(NTL::console::active).ntl;
		std::cout << "\n\nOptimising second NTL\n\n";
		ntl2 = opt2.optimise(NTL::console::active).ntl;

		std::cout << "d1:\t" << ntl1.get_d() << std::endl;
		std::cout << "d2:\t" << ntl2.get_d() << std::endl;
		NTL::NTL_sim sim;
		sim.add(ntl1);
		sim.add(ntl2);

		sim.w_h_profile(1, "W/H(z) (1)");
		sim.w_h_profile(2, "W/H(z) (2)");

		sim.set_f_sweep(1e7, 2.2e9);
		sim.s_matrix(1, 11, setup1.Zs, setup1.Zl, { "50", "100", "150" }, "S11 (1)");
		sim.s_matrix(2, 11, setup2.Zs, setup2.Zl, { "50", "25", "16.67" }, "S11 (2)");

		sim.merge("NTL");
		app.exec();
	}
	catch (std::invalid_argument& e)
	{
		std::cout << "Invalid argument exception: ";
		std::cout << e.what() << std::endl;
	}

	return 0;
}