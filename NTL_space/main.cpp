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
	QApplication app(argc, argv);

	NTL::NTL_opt_setup setup;
	setup.Z0 = 50; setup.Zs = 50; setup.Zl = 100; setup.er = 4.6; setup.d = 78e-3; setup.freqs = { 0.5e9, 1.3e9, 2e9 };
	setup.N = 10; setup.K = 50; setup.lb = std::vector<double>(setup.N, -1);  setup.ub = std::vector<double>(setup.N, 1);
	setup.toll_bounds = std::vector<double>(2, 1e-6); setup.toll_z = std::vector<double>(setup.K, 1e-6);
	setup.GBL_MAX = 350e3; setup.LCL_MAX = 50e3; setup.accepted_error = 6e-3;
	setup.m_Z_min = 0.1 * setup.Z0; setup.m_Z_max = 2.6 * setup.Z0; setup.max_attempts = 20;

	NTL::NTL_opt opt(setup);
	NTL::NTL ntl;

	try
	{
		ntl = opt.optimise(ntl, NTL::console::inactive);
	}
	catch (std::invalid_argument& e)
	{
		std::cout << e.what();
	}

	std::cout << "Result:\t" << ntl.get_Cn() << std::endl;

	NTL::NTL_sim sim(ntl);
	sim.set_f_sweep(1, 2e9);
	sim.z_profile();
	sim.w_h_profile();
	sim.s_matrix(setup.Zs, setup.Zl);

	sim.merge("NTL (1)");

	app.exec();

	return 0;
}